// 2D implicit FDM for steady state
// Neumann BC

#include <iostream>
#include "fdm.h"
#include <time.h>

extern void Gauss(double*,double*,double*,int);

int main()
{
  clock_t start, end;
  double cpuTime;
  std::ofstream aux_file;
  aux_file.open("cputime.txt");
  start = clock();
  //----------------------------------------------
  FDM* fdm = new FDM();
  fdm->SetInitialConditions();
  fdm->SetBoundaryConditions();
  //----------------------------------------------
  fdm->AssembleEquationSystem();
  Gauss(fdm->matrix,fdm->vecb,fdm->vecx,fdm->IJ);
  fdm->SaveTimeStep();
  fdm->OutputResults(0);
  //----------------------------------------------
  end = clock();
  cpuTime= (end-start)/ (double)(CLOCKS_PER_SEC);
  aux_file << "CPU time:" << cpuTime << std::endl;
  aux_file.close();
  fdm->out_file.close();
  fdm->eqs_file.close();
  delete fdm;
  return 0;
}
