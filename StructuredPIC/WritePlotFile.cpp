
#ifndef WRITEPLOTFILE_H
#define WRITEPLOTFILE_H

#include <string>
#include <vector>


#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_IntVect.H>

#include "WritePlotFile.h"


using namespace std;
using namespace amrex;

void WriteScalarData(const Vector<const MultiFab*>& rhs,
                     const Vector<const MultiFab*>& phi,
                     const Vector<Geometry>& geom, Real time, int step) //Vector<std::unique_ptr<EMParticleContainer>> particles)
{
    BL_PROFILE("WritePlotFile()");

    int num_levels = rhs.size();
    int num_output_comp = 2; 
    IntVect cc_flag = IntVect::TheZeroVector();
    Vector<unique_ptr<MultiFab>> output_cc(num_levels);
    int d_comp = 0; 
    for (int lev = 0; lev < num_levels; ++lev)
    {
        const BoxArray& nodal_ba = rhs[lev]->boxArray();
        output_cc[lev].reset(new MultiFab(convert(nodal_ba, cc_flag),
                                    rhs[lev]->DistributionMap(), num_output_comp, 0));
        average_node_to_cellcenter(*output_cc[lev], 0, *rhs[lev], 0, 1);
        average_node_to_cellcenter(*output_cc[lev], 1, *phi[lev], 0, 1);
    }

        Vector< string> varnames; 
        //scalar fields
        varnames.push_back("rhs");
        varnames.push_back("phi");

        Vector<int> level_steps;
        level_steps.push_back(0);
        level_steps.push_back(0);

        int output_levs = num_levels;

        Vector<IntVect> outputRR(output_levs);
        for (int lev = 0; lev < output_levs; ++lev) {
            outputRR[lev] = IntVect(D_DECL(2, 2, 2));
        }

        const std::string& pltfile = amrex::Concatenate("pltscalars", step, 5);
        WriteMultiLevelPlotfile(pltfile, output_levs, GetVecOfConstPtrs(output_cc),
                                varnames, geom, 0.0, level_steps, outputRR);
        Print() << "Done with plotfiles.";

};

void WriteVectorData(const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& BField,
                     const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& EField,
                     const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& jField,
                     const Vector<Geometry>& geom, Real time, int step)
{
    BL_PROFILE("WriteVectorData()");

    int num_levels = BField.size();
    int num_output_comp = 3*AMREX_SPACEDIM;
    
    IntVect cc_flag = IntVect::TheZeroVector(); 
    Vector<unique_ptr<MultiFab>> output_cc(num_levels);
    Vector<string> varnames;
    Vector<const MultiFab*> temp_fab(AMREX_SPACEDIM);

    int dcomp = 0; 
    for (int lev = 0; lev < num_levels; ++lev)
    {
         
        
        const BoxArray& cell_ba = BField[lev][0]->boxArray();
        output_cc[lev].reset(new MultiFab(convert(cell_ba, cc_flag),
                                    BField[lev][0]->DistributionMap(), num_output_comp, 0)); 
        
        temp_fab[0] = BField[lev][0];
        temp_fab[1] = BField[lev][1];
        temp_fab[2] = BField[lev][2];
        average_face_to_cellcenter(*output_cc[lev], dcomp, temp_fab);
        dcomp += 3; 
        
        temp_fab[0] = EField[lev][0];
        temp_fab[1] = EField[lev][1];
        temp_fab[2] = EField[lev][2];
        average_edge_to_cellcenter(*output_cc[lev], dcomp, temp_fab);
        dcomp += 3; 

        temp_fab[0] = jField[lev][0];
        temp_fab[1] = jField[lev][1];
        temp_fab[2] = jField[lev][2];
        average_edge_to_cellcenter(*output_cc[lev], dcomp, temp_fab);
        dcomp += 3; 

    }

    varnames.push_back("Bx");
    varnames.push_back("By");
    varnames.push_back("Bz");

    varnames.push_back("Ex");
    varnames.push_back("Ey");
    varnames.push_back("Ez");
    
    varnames.push_back("jx");
    varnames.push_back("jy");
    varnames.push_back("jz");

    int output_levs = num_levels;

    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputRR[lev] = IntVect(D_DECL(2, 2, 2));
    }

    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    const std::string& pltfile = amrex::Concatenate("pltvectors", step, 5);

    WriteMultiLevelPlotfile(pltfile, output_levs, GetVecOfConstPtrs(output_cc),
                                varnames, geom, 0.0, level_steps, outputRR);
    Print() << "Done with plotfiles.";
    
};

#endif
