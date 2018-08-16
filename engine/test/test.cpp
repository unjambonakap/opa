#include <gtest/gtest.h>

#include <opa_common.h>
#include <opa/engine/common.h>
#include <opa/engine/MipMapBuilder.h>
#include <opa/engine/gdal_obj.h>

using namespace std;
using namespace opa::math::common;
using namespace opa::engine;

TEST(ImageSampler, Test1) {
    Image img;
    int w = 16;
    int h = 16;
    img.init(w, h, Image::ImageFormat::Grey);
    REP(i, w) REP(j, w) img.set_col(i, j, Color(15 * (i / 2 + j / 2), 0, 0));
    ImageResampler sampler;
    OPA_DISP("orig", img);
    int factor=8;
    ImagePtr res = sampler.resample(
        ImageResampler::Params({ ImageResampler::ImageEntry(&img, 0, 0) }), w / factor, h / factor);
    OPA_DISP("img >> ", *res);
}

TEST(GdalSampler, Test2){

    GdalLoader loader;
    //loader.init(tiff_file);
    //Poco::File(dest).createDirectories();


    //auto rasterband = loader.ds()->GetRasterBand(1);
    //ImagePtr img(loader.get_img(rasterband->GetXSize(), rasterband->GetYSize()));
    //MipMapBuilder builder;
    //builder.init(img.get(), dest);
    //builder.build();

}

GTEST_API_ int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    opa::init::opa_init(argc, argv);
    (void)RUN_ALL_TESTS();
}
