#include <fstream>
#include <iomanip>

#include <Geometry.H>
#include <VisMF.H>

void writePlotFile (const std::string& dir, MultiFab& mf, std::vector<std::string> varnames,
Geometry& geom,Real& time);
