#ifndef FDM_H
#define FDM_H

#include <vector>
#include <list>
#include <fstream>
#include <iostream>

class BC
{
private:
public:
    BC();
    int bc_type;
    int geo_type;
    int node_number;
    double value;
};

class FDM
{
private:
public:
  int n,nn;
  int i,j;
  int I;
  int J;
  int IJ;
  int IJIJ;
  int i_row,j_col;
  double u0;
  double dx,dy;
  double x,y,x0,y0;
  double dx2;
  double dy2;
  double* matrix;
  double* vecb;
  double* vecx;
  std::vector<double>u;
  std::vector<int>bc_nodes;
  std::vector<double>u_bc;
  std::vector<BC*>bc_neumann;
  std::ofstream out_file;
  std::ofstream eqs_file;
  std::ofstream bc_file;
  BC* bc;
public:
    FDM();
    ~FDM();
    void SetInitialConditions();
    void SetBoundaryConditions();
    void SaveTimeStep();
    void OutputResults(int);
    void OutputMesh();
    void AssembleEquationSystem();
    void WriteEquationSystem();
    void IncorporateDirichletBoundaryConditions();
    void IncorporateNeumannBoundaryConditions();
    bool IsNodeNeumannBoundaryCondition(BC*);
};

#endif // FDM_H
