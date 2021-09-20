#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_WriteEBSurface.H>

#include <vector>

using namespace std; 
using namespace amrex; 

void torus(); 

int main(int argc, char* argv[]){

    amrex::Initialize(argc, argv);
    torus();
    amrex::Finalize();

}

void torus(){

    double R0 = 0.2477; 
    double a = 0.09; 
    int n_cell = 128; 
    int max_grid_size = 32; 

    Geometry geom;
    BoxArray grids;
    DistributionMapping dmap;
        RealBox rb({-0.35,-0.35,-0.15}, {0.35,0.35,0.15}); // physical domain
        Array<int,AMREX_SPACEDIM> is_periodic{false, false, false};
        Geometry::Setup(&rb, 0, is_periodic.data());
        Box domain(IntVect(0), IntVect(n_cell-1));
        geom.define(domain);
        grids.define(domain);
        grids.maxSize(max_grid_size);

        dmap.define(grids);

    EB2::TorusIF torus(Real(R0), Real(a), {0.0,0.0,0.0} ,true);
    auto gshop = EB2::makeShop(torus);
    EB2::Build(gshop,geom,0,0); 

    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > factory =
           makeEBFabFactory(geom, grids, dmap, {4, 4, 2}, EBSupport::full);
        const EBFArrayBoxFactory* ebfact = &(static_cast<amrex::EBFArrayBoxFactory const&>(*factory));

    MultiFab mf;
    {
        BoxArray ba(geom.Domain());
        ba.maxSize(max_grid_size);
        DistributionMapping dm{ba};

        std::unique_ptr<EBFArrayBoxFactory> factory
            = amrex::makeEBFabFactory(geom, ba, dm, {2,2,2}, EBSupport::volume);

        mf.define(ba, dm, 1, 0, MFInfo(), *factory);
        mf.setVal(1.0);
    }

    EB_WriteSingleLevelPlotfile("plt", mf, {"rho"}, geom, 0.0, 0);
    WriteEBSurface (grids, dmap, geom, ebfact);
}
