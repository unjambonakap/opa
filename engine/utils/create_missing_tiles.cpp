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

const float eps = 1e-3;

DEFINE_string(files, "", "");
DEFINE_string(dest_dir, "", "");

struct vec_cmp2 {

    bool operator()(const float &a, const float &b) const { return OPA_FLOAT_LT(a, b, eps); }
};
std::map<float, int, vec_cmp2> mpx;
std::map<float, int, vec_cmp2> mpy;

pii get(Pos2 a) { return pii(mpx[a.x], mpy[a.y]); }

void update_best(pii &best, pii &best2, pii a, pii b) {
    pii tmp = b - a;
    pii tmp2 = best2 - best;
    if (tmp.ST * tmp.ND > tmp2.ST * tmp2.ND)
        best = a, best2 = b;
}

int main(int argc, char **argv) {
    opa::init::opa_init(argc, argv);
    OPA_CHECK0(FLAGS_files.size() > 0);
    OPA_CHECK0(FLAGS_dest_dir.size() > 0);

    std::set<std::string> set_files;
    Poco::Glob::glob(FLAGS_files, set_files);
    std::vector<std::string> files(ALL(set_files));

    std::vector<proto::GdalData> data;
    for (auto &file : files) {
        OPA_CHECK0(Poco::File(file).exists());
        std::ifstream ifs(file, std::ifstream::binary);
        data.pb(proto::GdalData());
        OPA_CHECK0(data.back().ParseFromIstream(&ifs));
        data.back().clear_mmqt();
    }

    std::vector<Pos2> lst;
    for (auto &x : data) {
        Pos2 base = Pos2(x.x(), x.y());
        Pos2 dp = Pos2(x.dx() * (x.nx() + 1), x.dy() * (x.ny() + 1));
        for (auto p : UnitSquare)
            lst.pb(base + dp * p);
    }

    lst.pb(Pos2(-PI, -PI / 2));
    lst.pb(Pos2(PI, PI / 2));

    for (auto &u : lst) {
        mpx[u.x] = 0;
        mpy[u.y] = 0;
    }
    vector<float> rmpx(mpx.size());
    vector<float> rmpy(mpy.size());

    int nx = 0;
    int ny = 0;
    for (auto &a : mpx)
        rmpx[a.ND = nx++] = a.ST;
    for (auto &a : mpy)
        rmpy[a.ND = ny++] = a.ST;

    vector<vi> used;
    REP(i, nx) used.pb(vi(ny, -2));

    REP(i, data.size()) {
        vector<pii> tb;
        REP(j, 4) tb.pb(get(lst[4 * i + j]));
        pii lf = *min_element(ALL(tb));
        pii ur = *max_element(ALL(tb));
        FOR(x, lf.ST, ur.ST) FOR(y, lf.ND, ur.ND) used[x][y] = -1;
    }

    REP(i, nx) REP(j, ny) {
        if (used[i][j] == -1)
            continue;
        int oj = j;
        for (; used[i][j] == -2 && j < ny; ++j)
            ;

        for (; oj < j; ++oj)
            used[i][oj] = j;
        --j;
    }

    vector<pair<pii, pii> > created;
    REP(j, ny) REP(i, nx) {
        if (used[i][j] == -1)
            continue;
        vector<pii> tb; // start,height
        pii best = MP(i, j);
        pii best_end = best - MP(0, 0);
        tb.pb(MP(i, used[i][j]));

        for (;; ++i) {
            int h = i == nx ? -1 : used[i][j];
            while (tb.size() && tb.back().ND > h) {
                pii back = tb.back();
                tb.pop_back();
                update_best(best, best_end, MP(back.ST, j), MP(i, back.ND));
            }
            if (h == -1)
                break; // hack
        }
        created.pb(MP(best, best_end));
        FOR(x, best.ST, best_end.ST) FOR(y, best.ND, best_end.ND) {
            OPA_CHECK0(used[x][y] != -1);
            used[x][y] = -1;
        }
        i = best_end.ST - 1;
    }

    REP(i, nx) REP(j, ny) OPA_CHECK(used[i][j] == -1);
    Poco::File(FLAGS_dest_dir).createDirectories();

    ResourceManager rsc(ResourceManager::Params("created_", ".info", FLAGS_dest_dir));
    for (auto x : created) {
        proto::GdalData u;
        OPA_DISP("Create ", x);
        u.set_x(rmpx[x.ST.ST]);
        u.set_y(rmpy[x.ST.ND]);
        u.set_dx(rmpx[min(nx - 1, x.ND.ST)] - rmpx[x.ST.ST]);
        u.set_dy(rmpy[min(ny - 1, x.ND.ND)] - rmpy[x.ST.ND]);
        u.set_a1(data[0].a1());
        u.set_a2(data[0].a2());
        u.PrintDebugString();

        u.set_nx(1);
        u.set_ny(1);
        Mmqt<int> mmqt;
        Mmqt<int>::Node &node = mmqt.pool().get_new2();
        node.init();
        mmqt.root() = node.id;
        node.w = 1;
        node.h = 1;
        node.img_w = 1;
        node.img_h = 1;
        node.xo = 0;
        node.yo = 0;
        mmqt.store(*u.mutable_mmqt());

        std::ofstream ofs(rsc.create().ST.path(), std::ofstream::binary);
        u.SerializeToOstream(&ofs);
    }

    return 0;
}
