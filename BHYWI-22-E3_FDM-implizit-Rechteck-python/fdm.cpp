#include "fdm.h"
#include <stdlib.h>

enum {N,E,S,W,NE,SE,SW,NW};

FDM::FDM()
{
  I = 21; //43;
  J = 11; //33;
  IJ = I*J;
  IJIJ = IJ*IJ;
  dx = 10.;
  dy = 10.;
  u0 = 0.;
  x0 = 0.; //4423656.0991422;
  y0 = 0.; //5716944.1927754;
  //memory allocation
  u.resize(IJ);
  //output
  out_file.open("out.txt");
  eqs_file.open("matrices.txt");
  bc_file.open("bc.txt");
  dx2 = dx*dx;
  dy2 = dy*dy;
  matrix = new double[IJ*IJ];
  vecb = new double[IJ];
  vecx = new double[IJ];
}
FDM::~FDM()
{
  delete [] matrix;
  delete [] vecb;
  delete [] vecx;

  while(bc_neumann.size()>0)
  {
     delete bc_neumann[bc_neumann.size()-1];
     bc_neumann.pop_back();
  }
}
void FDM::SetInitialConditions()
{
  for(j=0;j<J;j++)
  {
    nn = j*I;
    for( i=0;i<I;i++)
    {
      n = nn+i;
      u[n] = u0;
    }
  }
}

void FDM::SetBoundaryConditions()
{
  int l;
  //------------------------------------
  //Dirichlet
  for(int i=0;i<I;i++)
  {
    bc_nodes.push_back(i); u[i] = 1;  // S
    bc_file << "S: " << i << " value= " << u[i] << std::endl;
  }
  for(int i=0;i<I;i++)
  {
    l = I*(J-1)+i;
    bc_nodes.push_back(l); u[l] = -1; // N
    bc_file << "N: " << l << " value= " << u[l] << std::endl;
  }
  for(int j=1;j<J-1;j++)
  {
    // E
    l = I*j;
    bc_nodes.push_back(l); u[l] = 0;
    bc_file << "E: " << l << " value= " << u[l] << std::endl;
  }
  for(int j=1;j<J-1;j++)
  {
    // W
    l = I*j+I-1;
    bc_nodes.push_back(l); u[l] = 0;
    bc_file << "W: " << l << " value= " << u[l] << std::endl;
  }
/*
  //------------------------------------
  //Neumann
  for(int j=1;j<J-1;j++)
  {
    // E
    l = I*j;
    bc = new BC();
    bc->node_number = l;
    bc->bc_type = 2;
    bc->geo_type = E;
    bc->value = 0.0;
    bc_neumann.push_back(bc);
    bc_file << "E: " << bc->node_number << " value= " << bc->value << std::endl;
  }
  for(int j=1;j<J-1;j++)
  {
    // W
    l = I*j+I-1;
    bc = new BC();
    bc->node_number = l;
    bc->bc_type = 2;
    bc->geo_type = W;
    bc->value = 0.0;
    bc_neumann.push_back(bc);
    bc_file << "W: " << bc->node_number << " value= " << bc->value << std::endl;
  }
  //.......................
*/
}

void FDM::AssembleEquationSystem()
{
  // Matrix entries
  for(i_row=0;i_row<IJ;i_row++)
  {
    vecb[i_row] = u[i_row];
    for(j_col=0;j_col<IJ;j_col++)
    {
      n = i_row*IJ+j_col;
      matrix[n] = 0.0;
      if(i_row==j_col)
      {
        matrix[n] = 4.;
      }
      if(abs((i_row-j_col))==1)
        matrix[i_row*IJ+j_col] = -1.;
      if(abs((i_row-j_col))==I)
        matrix[i_row*IJ+j_col] = -1.;
    }
  }
  // Incorporate Neumann boundary conditions
  IncorporateNeumannBoundaryConditions();
  // Incorporate Dirichlet boundary conditions
  IncorporateDirichletBoundaryConditions();
  // Matrix output
  WriteEquationSystem();
}

void FDM::IncorporateDirichletBoundaryConditions()
{
  size_t i_bc;
  int i_row, k;
  for(i_bc=0;i_bc<bc_nodes.size();i_bc++)
  {
     i_row = bc_nodes[i_bc];
     // Null off-diangonal entries of the related row and columns
     // Apply contribution to RHS by BC
     for(j=0;j<IJ;j++)
     {
        if(i_row == j)
           continue; // do not touch diagonals
        matrix[i_row*(IJ)+j] = 0.0; // NUll row
        k = j*(IJ)+i_row;
        // Apply contribution to RHS by BC
        vecb[j] -= matrix[k]*u[i_row];
        matrix[k] = 0.0; // Null column
     }
     // Apply Dirichlet BC
     vecb[i_row] = u[i_row]*matrix[i_row*(IJ)+i_row];
  }
}

void FDM::IncorporateNeumannBoundaryConditions()
{
  for(int i=0;i<(int)bc_neumann.size();i++)
  {
    bc = bc_neumann[i];
    n = bc->node_number;
    switch(bc->geo_type)
    {
      case S: //du/dy=0, South
        matrix[n*IJ+n+I] = -2.0;
        break;
      case W: //du/dx=0, West
        matrix[n*IJ+n+1] = -2.0;
        matrix[n*IJ+n-1] =  0.0;
        break;
      case N: //du/dy=0, North
        matrix[n*IJ+n-I] = -2.0;
        break;
      case E: //du/dx=0, East
        matrix[n*IJ+n-1] = -2.0;
        matrix[n*IJ+n+1] =  0.0;
        break;
      case SW: //du/dx=0,
        matrix[n*IJ+n+1] = -2.0;
        matrix[n*IJ+n+I] = -2.0;
        break;
      case NW: //du/dx=0,
        matrix[n*IJ+n+1] = -2.0;
        matrix[n*IJ+n-I] = -2.0;
        matrix[n*IJ+n-1] =  0.0;
        break;
      case NE: //du/dx=0
        matrix[n*IJ+n-1] = -2.0;
        matrix[n*IJ+n-I] = -2.0;
        break;
      case SE: //du/dx=0,
        matrix[n*IJ+n-1] = -2.0;
        matrix[n*IJ+n+I] = -2.0;
        matrix[n*IJ+n+1] =  0.0;
        break;
    }
  }
}

void FDM::SaveTimeStep()
{
  //save time step
  for(int j=0;j<J;j++)
    for(int i=0;i<I;i++)
    {
      u[j*I+i] = vecx[j*I+i];
    }
}

void FDM::OutputResults(int t)
{
  if((t%1)==0)
  {
    //out_file << "VARIABLES = X,Y,H"  << std::endl;
    out_file << "#ZONE T=\"BIG ZONE\", I=" << I << ", J=" << J <<", DATAPACKING=POINT"  << std::endl;
    for(int j=0;j<J;j++) //y
    {
      y = y0 + j*dy;
      nn = j*I;
      for(int i=0;i<I;i++) //x
      {
        n = nn+i;
        x = x0 + i*dx;
        out_file << x << "\t" << y << "\t" << u[n] << std::endl;
      }
    }
  }
}

void FDM::OutputMesh()
{
  //-----------------------------------------------------------------------
  //output:MSH
  std::ofstream msh_file;
  msh_file.open("2dfd.msh");
  msh_file << "#FEM_MSH" << std::endl;
  msh_file << " $PCS_TYPE" << std::endl;
  msh_file << "  GROUNDWATER_FLOW" << std::endl;
  msh_file << " $NODES" << std::endl;
  msh_file << "  " << IJ << std::endl;
  for(int j=0;j<J;j++)
  {
    y = y0 + j*dy;
    nn = j*I;
    for(int i=0;i<I;i++)
    {
      n = nn+i;
      x = x0 + i*dx;
      msh_file << n << "\t" << x << "\t" << y << "\t" << 0 << std::endl;
    }
  }
  msh_file << " $ELEMENTS" << std::endl;
  msh_file << "  " << (I-1)*(J-1) << std::endl;
  int c=0;
  for(int j=0;j<J-1;j++)
  {
    nn = j*I;
    for(int i=0;i<I-1;i++)
    {
      n = nn+i;
      msh_file << c++ << "\t" << 0 << " quad " << n << "\t" << n+1 << "\t" << n+I+1 << "\t" << n+I << std::endl;
    }
  }
  msh_file << "#STOP" << std::endl;
}

void FDM::WriteEquationSystem()
{
  for(i=0;i<IJ;i++)
  {
    for(j=0;j<IJ;j++)
    {
      eqs_file << matrix[i*IJ+j] << "\t";
    }
    eqs_file << "b:" << vecb[i] << " ";
    eqs_file << std::endl;
  }
}

bool FDM::IsNodeNeumannBoundaryCondition(BC* bc)
{
  bool is_node_bc = false;
  size_t k;          // site_t is the unsigned integer type, may not good for big size problems
  for(k=0;k<(size_t)bc_neumann.size();k++)
  {
    bc = bc_neumann[k];
    if(n==bc->node_number)
    {
      is_node_bc = true;
      return is_node_bc;
    }
  }
  return is_node_bc;
}

BC::BC()
{
}

