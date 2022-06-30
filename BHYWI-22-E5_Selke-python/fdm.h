#ifndef FDM_H
#define FDM_H

#include <vector>
#include <list>
#include <fstream>
#include <iostream>

class FDM
{
private:
public:
  //data structures
  std::vector<float>u_new;
  std::vector<float>u;
  std::vector<float>u_bc;
  float u0;
  float dx,dy,dt;
  float S0,Kf,Q;
  long i, j;
  int ix;
  int jy;
  std::ofstream out_file;
  std::fstream vtk_file;
  std::vector<int>bc_nodes;
  std::vector<int>nodes_inactive;
  int n,nn;
  float x,y,x0,y0;
  std::list<int>nodes_active;
  std::ifstream active_nodes_file;
  float dx2;
  float dy2;
  bool IsBCNode(int,std::vector<int>);
  bool IsNodeInactive(int,std::vector<int>);
  bool NodeInList(int,std::list<int>);
public:
    FDM();
    void SetActiveNodes();
    void SetInactiveNodes();
    void SetInitialConditions();
    void SetBoundaryConditions();
    void RunTimeStep();
    void SaveTimeStep();
    void OutputResults(int);
    void OutputMesh();
    void OutputResultsVTK(int);
};

#endif // FDM_H
