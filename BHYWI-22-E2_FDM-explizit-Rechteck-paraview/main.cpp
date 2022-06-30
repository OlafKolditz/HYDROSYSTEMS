#include <vector>
#include <fstream>
#include <time.h>

bool IsBCNode(int,std::vector<int>);
bool IsNodeInactive(int,std::vector<int>);

//int main(int argc, char *argv[])
int main()
{
  //clock_t start, end;
  double start, end;
  double cpuTime;
  //data structures
  std::vector<double>u_new;
  std::vector<double>u;
  std::vector<double>u_bc;
  double u0;
  double dx,dy,dt;
  double S0,Kf,Q;
  int ix;
  int jy;
  std::ofstream out_file;
  std::ofstream aux_file;
  std::vector<int>bc_nodes;
  std::vector<int>nodes_inactive;
  int n,nn;
  double x,y,x0,y0;
  //output
  out_file.open("out.dat");
  aux_file.open("cputime.txt");
  //-----------------------------------------------------------------------
  //parameters
  ix = 21;
  jy = 11;
  dx = 10.;
  dy = 10.;
  dt = 25.;//0.5e2; // sec
  //dt = 50.; // sec
  S0 = 1.e-5; // 1/m
  Kf = 1.e-5; // m/s
  Q = 0.0;
  //Q = 3.171e-9; // 100mm/year : 0.1/365/86400 = 3.171e-9 m/s
  u0 = 0.;
  x0 = 0.; //4423656.0991422;
  y0 = 0.; //5716944.1927754;
  //memory allocation
  u.resize(ix*jy);
  u_new.resize(ix*jy);
  //-----------------------------------------------------------------------
  //geometry
  //point_in_surface method
  //inactivate nodes
  //nodes_inactive.push_back(44);
  //nodes_inactive.push_back(1374);
  //-----------------------------------------------------------------------
  //initial conditions
  for(int i=0;i<ix;i++)
    for(int j=0;j<jy;j++)
    {
      u[j*(ix+1)] = u0;
      u_new[j*(ix+1)] = u0;
    }
  //-----------------------------------------------------------------------
  //boundary conditions
   // BC nodes from input file
  //top and bottom
  int l;
  for(int i=0;i<ix;i++)
  {
    bc_nodes.push_back(i); u[i] = 1.0;  u_new[i] = 1.0;
    l = ix*(jy-1)+i;
//aux_file << "l:" << l << std::endl;
    bc_nodes.push_back(l); u[l] = -1.0;  u_new[l] = -1.0;
  }
  //left and right side
  for(int j=1;j<jy-1;j++)
  {
    l = ix*j;
    bc_nodes.push_back(l); u[l] = 0.0;  u_new[l] = 0.0;
    l = ix*j+ix-1;
    bc_nodes.push_back(l); u[l] = 0.0;  u_new[l] = 0.0;
  }
  //optimization
  double dx2 = dx*dx;
  double dy2 = dy*dy;
  //=======================================================================
  int tn = 100;
  //output
  out_file << "TITLE = \"2D FDM\"" << std::endl;
  out_file << "VARIABLES = X, Y, H" << std::endl;
for(int t=0;t<tn;t++)
{
  start = clock();
  //-----------------------------------------------------------------------
  //time step
  for(int j=0;j<jy;j++)
  {
    nn = j*ix;
    for(int i=0;i<ix;i++)
    {
      n = nn+i;
      if(IsBCNode(n,bc_nodes))
        continue;
      //if(IsNodeInactive(n,nodes_inactive))
        //continue;
      u_new[n] = u[n] \
               + Kf/S0*dt/dx2 * (u[n+1]-2*u[n]+u[n-1]) \
               + Kf/S0*dt/dy2 * (u[(j+1)*ix+i]-2*u[n]+u[(j-1)*ix+i]) \
               + Q/S0;
    }
  }
  //-----------------------------------------------------------------------
  end = clock();
  // check steady state
  //-----------------------------------------------------------------------
  //save time step
  for(int j=0;j<jy;j++)
    for(int i=0;i<ix;i++)
    {
      u[j*ix+i] = u_new[j*ix+i];
    }
  //-----------------------------------------------------------------------
  //output:values
  out_file << "ZONE T=\"BIG ZONE\", I=" << ix << ", J=" << jy <<", DATAPACKING=POINT"  << std::endl;
  for(int j=0;j<jy;j++) //y
  {
    y = y0 + j*dy;
    nn = j*ix;
    for(int i=0;i<ix;i++) //x
    {
      n = nn+i;
      x = x0 + i*dx;
      if(IsNodeInactive(n,nodes_inactive))
        out_file << x << "\t" << y << "\t" << -9999 << std::endl;
      else
        out_file << x << "\t" << y << "\t" << u_new[n] << std::endl;
    }
  }
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
      msh_file << n << "\t" << x << "\t" << y << "\t" << 0 << std::endl;
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
      msh_file << c++ << "\t" << 0 << " quad " << n << "\t" << n+1 << "\t" << n+ix+1 << "\t" << n+ix << std::endl;
    }
  }
  msh_file << "#STOP" << std::endl;
}
  //=======================================================================
  //that's it
  cpuTime= (end-start)/ (double)(CLOCKS_PER_SEC);
  //cpuTime = difftime(end, start)/(double)CLOCKS_PER_SEC;
  aux_file << "CPU time:" << cpuTime << std::endl;
  aux_file.close();
  return 0;
}

bool IsBCNode(int n,std::vector<int>bc_nodes)
{
  bool is_node_bc = false;
  for(int k=0;k<(size_t)bc_nodes.size();k++)
  {
    if(n==bc_nodes[k])
    {
      is_node_bc = true;
      return is_node_bc;
    }
  }
  return is_node_bc;
}

bool IsNodeInactive(int n,std::vector<int>nodes_inactive)
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

  //std::ofstream xxx_file;
  //xxx_file.open("xxx.txt",std::ios::app);
//xxx_file << k << " " << n << " " << nodes_inactive[k] << std::endl;
