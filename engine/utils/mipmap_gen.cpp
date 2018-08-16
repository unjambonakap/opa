#include <opa_common.h>
#include <opa/engine/MipMapBuilder.h>
#include <opa/engine/gdal_obj.h>
#include <gdal.h>
#include <gdal_priv.h>
#include <Poco/Glob.h>
#include <Poco/Path.h>
#include <Poco/File.h>
#include <opa/utils/ProgressBar.h>

#include <opa/engine/proto/common.pb.h>
#include <opa/threading/runner.h>
#include <opa/threading/job.h>
#include <mipmap_gen.pb.h>
#include <opa/utils/ProgressBar.h>

using namespace std;
using namespace opa::engine;
using namespace opa::utils;
using namespace opa::threading;

DEFINE_string(tiff_example, "", "");
DEFINE_string(dest_dir, "", "");
DEFINE_bool(test_gdal, false, "");
DEFINE_bool(downscale, false, "");

void downscale(const std::string &file) {
    OPA_DISP0(file);
    OPA_CHECK0(Poco::File(file).exists());
    std::ifstream ifs(file, std::ifstream::binary);
    proto::GdalData gdal_proto;
    OPA_CHECK0(gdal_proto.ParseFromIstream(&ifs));
    Mmqt<int> mmqt;
    mmqt.load(gdal_proto.mmqt());

    int dest_side = 32;
    std::string path = Poco::Path(file).parent().toString();
    ResourceManager rsc(ResourceManager::Params(MM_RSC_PREFIX, "", path));

    while (1) {
        Mmqt<int>::Node &destnode = mmqt.pool().get_new2();
        Mmqt<int>::Node &node = mmqt.get(mmqt.root());
        if (node.img_w <= dest_side) {
            mmqt.pool().remove(destnode.id);
            break;
        }

        std::string filename = rsc.open(node.resource).path();
        Image img;
        std::ifstream tmp(filename, std::ifstream::binary);
        img.init(node.img_w, node.img_h, Image::ImageFormat::Grey);
        img.load(tmp);

        MmqtNodeId tmpid = destnode.id;
        destnode = node;
        destnode.id = tmpid;
        destnode.img_w /= 2;
        destnode.img_h /= 2;
        REP(i, 4) destnode.tb[i] = InvalidId;
        destnode.tb[0] = node.id;
        auto rscdata = rsc.create();
        destnode.resource = rscdata.ND;

        ImageResampler::Params params({ ImageResampler::ImageEntry(&img, 0, 0) });
        ImageResampler resampler;
        ImagePtr destimg = resampler.resample(params, destnode.img_w, destnode.img_h);
        std::ofstream ofs(rscdata.ST.path(), std::ofstream::binary);
        destimg->dump(ofs);

        destnode.tb[0] = node.id;
        mmqt.root() = destnode.id;
    }

    {
        mmqt.store(*gdal_proto.mutable_mmqt());
        std::ofstream ofs(file, std::ofstream::binary);
        gdal_proto.SerializeToOstream(&ofs);
    }
}

void build1(const std::string &tiff_file, const std::string &dest) {
    GdalLoader loader;
    loader.init(tiff_file);
    try {
        Poco::File(dest).createDirectories();
    } catch (...) {
    }
    Poco::File dest_file =
        Poco::File(Poco::Path(dest, Poco::Path(tiff_file).getFileName()).setExtension("info"));
    if (dest_file.exists() && dest_file.isFile())
        return;

    auto rasterband = loader.ds()->GetRasterBand(1);
    ImagePtr img = loader.get_img(rasterband->GetXSize(), rasterband->GetYSize());
    MipMapBuilder builder;
    builder.init(img.get(), dest);

    builder.build();

    proto::GdalData data;
    builder.get_proto(*data.mutable_mmqt());
    loader.setup_proto(data);

    std::string info_file = dest_file.path();
    std::ofstream os(info_file, std::ofstream::binary);
    data.SerializeToOstream(&os);
}

class GenJob : public TJob<mipmapgen::Empty, mipmapgen::Data, mipmapgen::Empty> {

  public:
    GenJob() {}
    GenJob(const vector<string> &files, const string &dest) {
        m_files = files;
        m_dest = dest;
    }
    void tworker_initialize(const JobMsg &base_msg, const mipmapgen::Empty &init) {}
    void tworker_do_work(const JobMsg &base_msg, const mipmapgen::Data &data,
                         mipmapgen::Empty &out_res) {
        if (data.dest().empty())
            downscale(data.filename());
        else
            build1(data.filename(), data.dest());
    }

    void tserver_initialize(mipmapgen::Empty &out_init) {}
    void
    tserver_get_work(const std::function<DataId(const mipmapgen::Data &data, bool &out_more)> &cb) {
        m_pbj = ProgressBar::get()->create(ProgressBar::JobDesc(m_files.size(), "mm gen"));
        mipmapgen::Data data;
        for (auto x : m_files) {
            data.set_filename(x);
            if (!m_dest.empty()) {
                auto filename = Poco::Path(x).setExtension("").getFileName();
                data.set_dest(Poco::Path(FLAGS_dest_dir, filename).toString());
            }

            bool more;
            cb(data, more);
        }
    }
    void tserver_set_work_result(const JobMsg &base_msg, const mipmapgen::Empty &res,
                                 DataId data_id) {
        m_pbj.increment(m_files[data_id]);
        ProgressBar::get()->refresh();
    }

  private:
    ProgressBar::ProgressJob m_pbj;
    vector<string> m_files;
    string m_dest;
};

int main(int argc, char **argv) {
    opa::init::opa_init(argc, argv);

    if (FLAGS_test_gdal) {

        if (FLAGS_downscale) {
            downscale(FLAGS_tiff_example);

        } else {

            OPA_CHECK0(FLAGS_tiff_example.size() != 0);
            OPA_CHECK0(FLAGS_dest_dir.size() != 0);
            GdalLoader loader;
            loader.init(FLAGS_tiff_example);

            auto rasterband = loader.ds()->GetRasterBand(1);
            ImagePtr img(loader.get_img(rasterband->GetXSize(), rasterband->GetYSize()));

            ImageResampler sampler;
            int factor = 6;
            ImagePtr res = sampler.resample(
                ImageResampler::Params({ ImageResampler::ImageEntry(img.get(), 0, 0) }),
                img->w() / factor, img->h() / factor);

            std::ofstream ofs(FLAGS_dest_dir, std::ofstream::binary);
            res->dump(ofs);
        }
    } else {

        if (FLAGS_downscale) {

            std::set<std::string> set_files;
            Poco::Glob::glob(FLAGS_tiff_example, set_files);
            std::vector<std::string> files(ALL(set_files));

            string name = "1";
            auto id = Runner::Register_job(name, []() { return new GenJob; });
            Runner runner;
            runner.run_both();
            Dispatcher *x = runner.dispatcher();
            GenJob job(files, "");
            puts("go process");
            x->process_job(job, id);
            puts("DONE JOB");
        } else {

            std::set<std::string> set_files;
            Poco::Glob::glob(FLAGS_tiff_example, set_files);
            std::vector<std::string> files(ALL(set_files));

            if (0) {

                ProgressBar *pb = ProgressBar::get();
                auto job = pb->create(ProgressBar::JobDesc(files.size(), "gdal mipmap builder"));

                for (auto &x : files) {
                    auto filename = Poco::Path(x).setExtension("").getFileName();
                    build1(x, Poco::Path(FLAGS_dest_dir, filename).toString());

                    job.increment(format("done %s", filename.c_str()));
                    pb->refresh();
                }
            } else {

                string name = "1";
                auto id = Runner::Register_job(name, []() { return new GenJob; });
                Runner runner;
                runner.run_both();
                Dispatcher *x = runner.dispatcher();
                GenJob job(files, FLAGS_dest_dir);
                puts("go process");
                x->process_job(job, id);
                puts("DONE JOB");
            }
        }
    }

    return 0;
}
