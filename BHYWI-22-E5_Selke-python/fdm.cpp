#include "fdm.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

FDM::FDM()
{
  ix = 43;
  jy = 33;
  dx = 1000.;
  dy = 1000.;
  dt = 1e6; // 0.25e6 sec
  S0 = 1e-5;
  Kf = 1e-5;                                     // m/s
  //Q = 3.171e-9; // 100mm/year : 0.1/365/86400 = 3.171e-9 m/s
  Q = 0.; 
  u0 = 300.;
  x0 = 4423656.0991422;
  y0 = 5716944.1927754;
  //memory allocation
  u.resize(ix*jy);
  u_new.resize(ix*jy);
  //output
  out_file.open("out.txt");
  dx2 = dx*dx;
  dy2 = dy*dy;
}

bool FDM::IsBCNode(int n,std::vector<int>bc_nodes)
{
   bool is_node_bc = false;
   size_t k;          // site_t is the unsigned integer type, may not good for big size problems
   for(k=0;k<(size_t)bc_nodes.size();k++)
   {
      if(n==bc_nodes[k])
      {
         is_node_bc = true;
         return is_node_bc;
      }
   }
   return is_node_bc;
}

bool FDM::IsNodeInactive(int n,std::vector<int>nodes_inactive)
{
  bool is_node_inactive = false;
  for(int k=0;k<(size_t)nodes_inactive.size();k++)
  {
    if(n==nodes_inactive[k])
    {
      is_node_inactive = true;
      return is_node_inactive;
    }
  }
  return is_node_inactive;
}

bool FDM::NodeInList(int n,std::list<int>nodes_active)
{
  std::list<int>::const_iterator p = nodes_active.begin();
  while(p!=nodes_active.end())
  {
   if(n==*p)
     return true;
    ++p;
  }
  return false;
}

void FDM::SetActiveNodes()
{
  active_nodes_file.open("ActiveNodes.txt");
  if(!active_nodes_file.good())
  {
    std::cout << "Error: File not found" << std::endl;
    return;
  }
  int na;
  while(!active_nodes_file.eof())
  {
    active_nodes_file >> na;
    nodes_active.push_back(na);
  }
  nodes_active.sort();
  std::ofstream active_nodes_file_test;
  active_nodes_file_test.open("../ActiveNodesTest.txt");
  std::list<int>::const_iterator p = nodes_active.begin();
  while(p!=nodes_active.end())
  {
    active_nodes_file_test << *p << std::endl;
    ++p;
  }
  active_nodes_file_test.close();
}

void FDM::SetInactiveNodes()
{
  std::ofstream inactive_nodes_file;
  inactive_nodes_file.open("../InactiveNodes.txt");
  for(j=0;j<jy;j++)
  {
    nn = j*ix;
    for( i=0;i<ix;i++)
    {
      n = nn+i;
      if(!NodeInList(n,nodes_active))
        nodes_inactive.push_back(n);
    }
  }
  for(i=0;i<nodes_inactive.size();i++)
  {
    inactive_nodes_file << nodes_inactive[i] << std::endl;
  }
  inactive_nodes_file.close();
}

void FDM::SetInitialConditions()
{
  for(j=0;j<jy;j++)
  {
    nn = j*ix;
    for( i=0;i<ix;i++)
    {
      n = nn+i;
      u[n] = u0;
      u_new[n] = u0;
    }
  }
}

void FDM::SetBoundaryConditions()
{
  // BC nodes from input file
  //top and bottom
  int l;
  for( i=0;i<ix;i++)
  {
    bc_nodes.push_back(i); u[i] = u0;  u_new[i] = u0;
    l = ix*(jy-1)+i;
    if(l>1402&&l<1408)
    {
      bc_nodes.push_back(l); u[l] = u0;  u_new[l] = u0;
    }
    else
    {
      bc_nodes.push_back(l); u[l] = u0;  u_new[l] = u0;
    }
  }
  //left and right side
  for(j=1;j<jy-1;j++)
  {
    l = ix*j;
    if(j>4&&j<9)
    {
      bc_nodes.push_back(l); u[l] = 800.;  u_new[l] = 800.;
    }
    else
    {
      bc_nodes.push_back(l); u[l] = u0;  u_new[l] = u0;
    }
    l = ix*j+ix-1;
    bc_nodes.push_back(l); u[l] = u0;  u_new[l] = u0;
  }
}

void FDM::RunTimeStep()
{
  for( j=0;j<jy;j++)
  {
    nn = j*ix;
    for( i=0;i<ix;i++)
    {
      n = nn+i;
      if(IsBCNode(n,bc_nodes))
        continue;
      if(IsNodeInactive(n,nodes_inactive))
        continue;
      u_new[n] = u[n] \
               + Kf/S0*dt/dx2 * (u[n+1]-2*u[n]+u[n-1]) \
               + Kf/S0*dt/dy2 * (u[(j+1)*ix+i]-2*u[n]+u[(j-1)*ix+i]) \
               + Q/S0;
    }
  }
}

void FDM::SaveTimeStep()
{
  //save time step
  for(int j=0;j<jy;j++)
    for(int i=0;i<ix;i++)
    {
      u[j*ix+i] = u_new[j*ix+i];
    }
}

void FDM::OutputResults(int t)
{
//  if((t%10)==0)
  {
    out_file << "#ZONE T=\"BIG ZONE\", I=" << ix << ", J=" << jy <<", DATAPACKING=POINT"  << std::endl;
    for(int j=0;j<jy;j++) //y
    {
      y = y0 + j*dy;
      nn = j*ix;
      for(int i=0;i<ix;i++) //x
      {
        n = nn+i;
        x = x0 + i*dx;
        if(IsNodeInactive(n,nodes_inactive))
          out_file << x << "\t" << y << "\t" << 0.0 << std::endl;
        else
          out_file << x << "\t" << y << "\t" << u_new[n] << std::endl;
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
  msh_file << "  " << ix*jy << std::endl;
  for(int j=0;j<jy;j++)
  {
    y = y0 + j*dy;
    nn = j*ix;
    for(int i=0;i<ix;i++)
    {
      n = nn+i;
      x = x0 + i*dx;
      msh_file << n << " " << x << " " << y << " " << 0.0 << std::endl;
    }
  }
  msh_file << " $ELEMENTS" << std::endl;
  msh_file << "  " << (ix-1)*(jy-1) << std::endl;
  int c=0;
  for(int j=0;j<jy-1;j++)
  {
    nn = j*ix;
    for(int i=0;i<ix-1;i++)
    {
      n = nn+i;
      msh_file << c++ << " " << 0 << " quad " << n << " " << n+1 << " " << n+ix+1 << " " << n+ix << std::endl;
    }
  }
  msh_file << "#STOP" << std::endl;
}

void FDM::OutputResultsVTK(int t)
{
  std::string vtk_file_name("file");
  std::stringstream ss;
  ss << t;
  vtk_file_name += ss.str();
  vtk_file_name += ".vtk";
  std::fstream vtk_file (vtk_file_name.data(),std::ios::out);
  vtk_file.setf(std::ios::scientific,std::ios::floatfield);
  vtk_file.precision(12);
  if (!vtk_file.good())
  {
    std::cout << "Could not open file for writing: VTK " << std::endl;
    return;
  }
  vtk_file.seekg(0L,std::ios::beg);
  //Header
  vtk_file << "# vtk DataFile Version 3.0" << std::endl;
  vtk_file << "Unstructured Grid from GW3" << std::endl;
  vtk_file << "ASCII" << std::endl;
  vtk_file << "DATASET UNSTRUCTURED_GRID" << std::endl;
  //Point coordinates
  vtk_file << "POINTS " << ix*jy << " double" << std::endl;
  for(int j=0;j<jy;j++)
  {
    y = y0 + j*dy;
    nn = j*ix;
    for(int i=0;i<ix;i++)
    {
      n = nn+i;
      x = x0 + i*dx;
      vtk_file << x << " " << y << " " << 0.0 << std::endl;
    }
  }
  //
  // count overall length of element vector
  long numAllPoints =0;
  for(size_t i=0; i < (ix-1)*(jy-1); i++)
  {
    numAllPoints = numAllPoints + 4 + 1;
  }

  //Cells
  vtk_file << "CELLS " << (ix-1)*(jy-1) << " " << numAllPoints << std::endl;
  for(int j=0;j<jy-1;j++)
  {
    nn = j*ix;
    for(int i=0;i<ix-1;i++)
    {
      n = nn+i;
      vtk_file << "4 " << n << " " << n+1 << " " << n+ix+1 << " " << n+ix << std::endl;
    }
  }
  vtk_file << std::endl;
  //Cell types
  vtk_file << "CELL_TYPES " << (ix-1)*(jy-1) << std::endl;
  for(int j=0;j<jy-1;j++)
  {
    for(int i=0;i<ix-1;i++)
    {
      vtk_file << "9 " << std::endl;
    }
  }
  vtk_file << std::endl;
  //Point data
  vtk_file << "POINT_DATA " << ix*jy << std::endl;
  vtk_file << "SCALARS HEAD double 1" << std::endl;
  vtk_file << "LOOKUP_TABLE default" << std::endl;
  for(int j=0;j<jy;j++)
  {
    y = y0 + j*dy;
    nn = j*ix;
    for(int i=0;i<ix;i++)
    {
      n = nn+i;
      vtk_file << u_new[n] << std::endl;
    }
  }
}


