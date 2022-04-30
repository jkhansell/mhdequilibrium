
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Print.H>

#ifdef AMREX_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#include "AmrCoreAdv.H"
#include "Prob.H"
#include "solctra_utils.H"
#include "Tagging.H"

using namespace amrex;

// constructor - reads in parameters from inputs file
//             - sizes multilevel arrays and data structures
//             - initializes BCRec boundary condition object
AmrCoreAdv::AmrCoreAdv(GlobalData& coils, Coil* Rmi, Coil* Rmf, Real r, Real deltah)
{
    ReadParameters();

    // Geometry on all levels has been defined already.

    // No valid BoxArray and DistributionMapping have been defined.
    // But the arrays for them have been resized.

    Coils = coils; 
    rmi = Rmi; 
    rmf = Rmf;
    R0 = r; 
    dh = deltah;

    int nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    /*
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }*/

    t.resize(nlevs_max, 0.0);
    dt.resize(nlevs_max, 1.e100);

    Bx.resize(nlevs_max); 
    By.resize(nlevs_max); 
    Bz.resize(nlevs_max); 

    gradBx.resize(nlevs_max); 
    gradBy.resize(nlevs_max); 
    gradBz.resize(nlevs_max); 

    Jacobian.resize(nlevs_max); 

    // periodic boundaries
    int bc_lo[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};
    int bc_hi[] = {BCType::int_dir, BCType::int_dir, BCType::int_dir};

    bcs.resize(1);     // Setup 1-component
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        // lo-side BCs
        if (bc_lo[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_lo[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_lo[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setLo(idim, bc_lo[idim]);
        }
        else {
            amrex::Abort("Invalid bc_lo");
        }

        // hi-side BCSs
        if (bc_hi[idim] == BCType::int_dir  ||  // periodic uses "internal Dirichlet"
            bc_hi[idim] == BCType::foextrap ||  // first-order extrapolation
            bc_hi[idim] == BCType::ext_dir ) {  // external Dirichlet
            bcs[0].setHi(idim, bc_hi[idim]);
        }
        else {
            amrex::Abort("Invalid bc_hi");
        }
    }


}

AmrCoreAdv::~AmrCoreAdv ()
{
}

// initializes multilevel data
void
AmrCoreAdv::InitData ()
{
    const Real time = 0.0;
    InitFromScratch(time);
    AverageDown();
    
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    const int ncomp = Bx[lev-1].nComp();
    const int ncompjacobian = Jacobian[lev-1].nComp();
    const int nghost = Bx[lev-1].nGrow();
    const int nghostjacobian = Bx[lev-1].nGrow();


    Bx[lev].define(ba, dm, ncomp, nghost);
    By[lev].define(ba, dm, ncomp, nghost);
    Bz[lev].define(ba, dm, ncomp, nghost);
    gradBx[lev].define(ba, dm, ncomp, nghost);
    gradBy[lev].define(ba, dm, ncomp, nghost);
    gradBz[lev].define(ba, dm, ncomp, nghost);

    Jacobian[lev].define(ba, dm, ncompjacobian, nghostjacobian);

    t[lev] = time;

    FillCoarsePatch(lev, time, Bx[lev], 0, ncomp);
    FillCoarsePatch(lev, time, By[lev], 0, ncomp);
    FillCoarsePatch(lev, time, Bz[lev], 0, ncomp);
    FillCoarsePatch(lev, time, gradBx[lev], 0, ncomp);
    FillCoarsePatch(lev, time, gradBy[lev], 0, ncomp);
    FillCoarsePatch(lev, time, gradBz[lev], 0, ncomp);
    FillCoarsePatch(lev, time, Jacobian[lev], 0, ncompjacobian);
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    const int ncomp = Bx[lev-1].nComp();
    const int ncompjacobian = Jacobian[lev-1].nComp();
    const int nghost = Bx[lev-1].nGrow();
    const int nghostjacobian = Bx[lev-1].nGrow();

    MultiFab new_state_field(ba, dm, ncomp, nghost);

    MultiFab new_state_jacobian(ba, dm, ncomp, nghostjacobian);

    FillPatch(lev, time, new_state_field, 0, ncomp);
    FillPatch(lev, time, new_state_field, 0, ncompjacobian);


    std::swap(new_state_field, Bx[lev]);
    std::swap(new_state_field, By[lev]);
    std::swap(new_state_field, Bz[lev]);

    std::swap(new_state_field, gradBx[lev]);
    std::swap(new_state_field, gradBy[lev]);
    std::swap(new_state_field, gradBz[lev]);
    
    std::swap(new_state_jacobian, Jacobian[lev]);
    
    t[lev] = time;

    // This clears the old MultiFab and allocates the new one
}

// Delete level data
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ClearLevel (int lev)
{
    Bx[lev].clear();
    By[lev].clear();
    Bz[lev].clear();

    gradBx[lev].clear();
    gradBy[lev].clear();
    gradBz[lev].clear();

    Jacobian[lev].clear();
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm)
{
    const int ncomp = 1;
    const int ncompjacobian = 9;
    const int nghost = 1;
    const int nghostjacobian = 1; 

    Bx[lev].define(ba, dm, ncomp, nghost);
    By[lev].define(ba, dm, ncomp, nghost);
    Bz[lev].define(ba, dm, ncomp, nghost);

    gradBx[lev].define(ba, dm, ncomp, nghost);
    gradBy[lev].define(ba, dm, ncomp, nghost);
    gradBz[lev].define(ba, dm, ncomp, nghost);

    Jacobian[lev].define(ba, dm, ncompjacobian, nghost);

    t[lev] = time;

    // This clears the old MultiFab and allocates the new one

    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();
    //Print() << problo[0] << problo[1];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(Bx[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();

        Array4<Real> iBx = Bx[lev][mfi].array();
        Array4<Real> iBy = By[lev][mfi].array();
        Array4<Real> iBz = Bz[lev][mfi].array();
        
        ParallelFor(box, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            init_fields_zero(i,j,k, iBx, iBy, iBz);
        });
    }
/*
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    for (MFIter mfi(Bx[lev], TilingIfNotGPU());mfi.isValid(); ++mfi) 
    {
        const Box& box = mfi.tilebox();
        
        Array4<Real> iBx = Bx[lev][mfi].array();
        Array4<Real> iBy = By[lev][mfi].array();
        Array4<Real> iBz = Bz[lev][mfi].array();

        Array4<Real> igradBx = gradBx[lev][mfi].array();
        Array4<Real> igradBy = gradBy[lev][mfi].array();
        Array4<Real> igradBz = gradBz[lev][mfi].array();

        Array4<Real> ijacobian = Jacobian[lev][mfi].array();

        ParallelFor(box, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            derivative_calculation(i, j, k, iBx, iBy, iBz, 
                                igradBx, igradBy, igradBz,
                                ijacobian, problo, dx);
        });
    }
*/
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    //const int clearval = TagBox::CLEAR;
    const int tagval = TagBox::SET;
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if(Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(Jacobian[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx  = mfi.tilebox();
            const auto tagfab  = tags.array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                state_error(i, j, k, tagfab, tagval, R0, dh, problo, dx);
                //Print() << i;
            });
        }
    }
}

// read in some parameters from inputs file
void
AmrCoreAdv::ReadParameters ()
{
    {
        ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }

    {
        ParmParse pp("amr"); // Traditionally, these have prefix, amr.

        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
    }
}

// set covered coarse cells to be the average of overlying fine cells
void
AmrCoreAdv::AverageDown ()
{
    for (int lev = finest_level-1; lev >= 1; --lev)
    {
        amrex::average_down(Bx[lev+1], Bx[lev],
                            geom[lev+1], geom[lev],
                            0, Bx[lev].nComp(), refRatio(lev));
        amrex::average_down(By[lev+1], By[lev],
                            geom[lev+1], geom[lev],
                            0, By[lev].nComp(), refRatio(lev));
        amrex::average_down(Bz[lev+1], Bz[lev],
                            geom[lev+1], geom[lev],
                            0, Bz[lev].nComp(), refRatio(lev));

        amrex::average_down(gradBx[lev+1], gradBx[lev],
                            geom[lev+1], geom[lev],
                            0, gradBx[lev].nComp(), refRatio(lev));
        amrex::average_down(gradBy[lev+1], gradBy[lev],
                            geom[lev+1], geom[lev],
                            0, gradBy[lev].nComp(), refRatio(lev));
        amrex::average_down(Bz[lev+1], Bz[lev],
                            geom[lev+1], geom[lev],
                            0, Bz[lev].nComp(), refRatio(lev));

        amrex::average_down(Jacobian[lev+1], Jacobian[lev],
                            geom[lev+1], geom[lev],
                            0, Jacobian[lev].nComp(), refRatio(lev));
    }
}

// more flexible version of AverageDown() that lets you average down across multiple levels
void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
        amrex::average_down(Bx[crse_lev+1], Bx[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, Bx[crse_lev].nComp(), refRatio(crse_lev));
        amrex::average_down(By[crse_lev+1], By[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, By[crse_lev].nComp(), refRatio(crse_lev));
        amrex::average_down(Bz[crse_lev+1], Bz[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, Bz[crse_lev].nComp(), refRatio(crse_lev));

        amrex::average_down(gradBx[crse_lev+1], gradBx[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, gradBx[crse_lev].nComp(), refRatio(crse_lev));
        amrex::average_down(gradBy[crse_lev+1], gradBy[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, gradBy[crse_lev].nComp(), refRatio(crse_lev));
        amrex::average_down(gradBz[crse_lev+1], Bz[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, Bz[crse_lev].nComp(), refRatio(crse_lev));

        amrex::average_down(Jacobian[crse_lev+1], Jacobian[crse_lev],
                            geom[crse_lev+1], geom[crse_lev],
                            0, Jacobian[crse_lev].nComp(), refRatio(crse_lev));
}

// compute a new multifab by coping in phi from valid region and filling ghost cells
// works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    if (lev == 0)
    {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bndry_func);
            amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
                                        geom[lev], physbc, 0);
        }
    }
    else
    {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        Interpolater* mapper = &cell_cons_interp;

        {
            CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
            PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
            PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

            amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
                                      0, icomp, ncomp, geom[lev-1], geom[lev],
                                      cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                      mapper, bcs, 0);
        }
    }
}

// fill an entire multifab by interpolating from the coarser level
// this comes into play when a new level of refinement appears
void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
    BL_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    Interpolater* mapper = &cell_cons_interp;

    if (cmf.size() != 1) {
        amrex::Abort("FillCoarsePatch: how did this happen?");
    }
    else
    {
        CpuBndryFuncFab bndry_func(nullptr);  // Without EXT_DIR, we can pass a nullptr.
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bndry_func);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bndry_func);

        amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
                                     cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                                     mapper, bcs, 0);
    }
}

// utility to copy in data from phi_old and/or phi_new into another multifab
void
AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    data.push_back(&Bx[lev]);
    data.push_back(&By[lev]);
    data.push_back(&Bz[lev]);
    data.push_back(&gradBx[lev]);
    data.push_back(&gradBy[lev]);
    data.push_back(&gradBz[lev]);
    data.push_back(&Jacobian[lev]);
 
    datatime.push_back(t[lev]);
}


// get plotfile name
std::string
AmrCoreAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

// put together an array of multifabs for writing
Vector<const MultiFab*>
AmrCoreAdv::PlotFileMF () const
{
    Vector<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
        r.push_back(&Bx[i]);
        r.push_back(&By[i]);
        r.push_back(&Bz[i]);
        r.push_back(&gradBx[i]);
        r.push_back(&gradBy[i]);
        r.push_back(&gradBz[i]);
        r.push_back(&Jacobian[i]);
    }
    return r;
}

// set plotfile variable names
Vector<std::string>
AmrCoreAdv::PlotFileVarNames () const
{
    return {"Bx", "By", "Bz", "gradBx", "gradBy", "gradBz", "Jacobian"};
}

// write plotfile to disk
void
AmrCoreAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
                                   Geom(), t[0], istep, refRatio());
}

