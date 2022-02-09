#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

#include "WritePlotFile.h"

using namespace amrex;

void WritePlotFile
(
    const Array<const MultiFab*, AMREX_SPACEDIM>& B,
    const Array<const MultiFab*, AMREX_SPACEDIM>& gradB,
    const MultiFab& jacobian, const Geometry& geom, 
    Real time, int step
)
{
    BL_PROFILE("WritePlotFile()"); 
    int num_output_comp = 6; 

    IntVect cc_flag = IntVect::TheZeroVector(); 
    
    std::unique_ptr<MultiFab> output_cc;
    Vector<std::string> varnames; 

    Vector<const MultiFab*> tempfab;
        
    int d_comp = 0; 

    const BoxArray& nodal_ba = B[0]->boxArray();
    output_cc.reset(new MultiFab(convert(nodal_ba, cc_flag), 
                        B[0]->DistributionMap(), num_output_comp, 0));

    average_node_to_cellcenter(*output_cc, d_comp, *B[0], 0, 1, 0);
    d_comp += 1;
    average_node_to_cellcenter(*output_cc, d_comp, *B[1], 0, 1, 0);
    d_comp += 1;
    average_node_to_cellcenter(*output_cc, d_comp, *B[2], 0, 1, 0);
    d_comp += 1;
    average_node_to_cellcenter(*output_cc, d_comp, *gradB[0], 0, 1, 0);
    d_comp += 1;
    average_node_to_cellcenter(*output_cc, d_comp, *gradB[1], 0, 1, 0);
    d_comp += 1;
    average_node_to_cellcenter(*output_cc, d_comp, *gradB[2], 0, 1, 0);
    d_comp += 1;

    varnames.push_back("Bx");
    varnames.push_back("By");
    varnames.push_back("Bz");
    varnames.push_back("gradBx");
    varnames.push_back("gradBy");
    varnames.push_back("gradBz");
    const std::string& pltfile = amrex::Concatenate("plotfiles/pltfields", step, 5);

    WriteSingleLevelPlotfile(pltfile, *output_cc,
                     varnames, geom, time, step);
}