#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>

//extern void TestOutput(double*);

int main(int argc, char *argv[])
{
  //-----------------------------------------------------------
  // 1 Datenstrukturen
  int n=6;
  int solver_iterations = 100;
  double eps = 1e-3;
  double x[n], x0[n];
  double y;
  for(int i=0;i<n;i++)
    x[i]=10000.;
  ofstream out_file;
  out_file.open("out.txt");
  //-----------------------------------------------------------
  // 2 Parameter
  double dx1,dx2,dx3,dx4,dx5,dx6,dy1,dy2,dy3,dy4,dy5,dy6;
  dx1 = 100.;
  dx2 = 1000.;
  dx3 = 1000.;
  dx4 = 100.;
  dx5 = 1000.;
  dx6 = 1000.;
  dy1 = 500.;
  dy2 = 500.;
  dy3 = 500.;
  dy4 = 500.;
  dy5 = 500.;
  dy6 = 500.;
  double T1,T2,T3,T4,T5,T6;
  double T12,T32,T52;
  T1 = 0.05;
  T2 = 0.05;
  T3 = 0.05;
  T4 = 0.01;
  T5 = 0.01;
  T6 = 0.01;
  T12 = (dx1+dx2)/(dx1/T1+dx2/T2);
  cout << "T12: " << T12 << endl;
  T32 = (dx3+dx2)/(dx3/T3+dx2/T2);
  cout << "T32: " << T32 << endl;
  T52 = (dy5+dy2)/(dy5/T5+dy2/T2);
  cout << "T52: " << T52 << endl;
  //Randbedingungen
  double h1,h4;
  h1 = 15.;
  h4 = 15.;
  double Q12,Q32,Q52,QR,QP2,QP3,QP5,QP6;
  // QR = R*dx2*dy2
  QR = 1.e-07 * dx2 * dy2;
  cout << "QR: " << QR << endl;
  QP2 = - 0.1;
  QP3 = 0.0;
  QP5 = 0.0;
  //  QP6 = - 0.005;
  QP6 = -0.1;
  //-----------------------------------------------------------
  // 2 Knotenbilanzen
  //...........................................................
  // 2.2: Q12+Q32+Q52+QR+QP2=0
  cout << "Knoten 2:" << endl;
  //...........................................................
  // Q12 = dy12 * T12 *(h1-h2)/(dx1/2+dx2/2)
  // Q12 = c121*h1 + c122*h2
  double c121, c122;
  c121 = dy1 * T12 * 2./(dx1+dx2);
  cout << "c121: " << c121 << endl;
  c122 = - c121;
  cout << "c122: " << c122 << endl;
  // Q32 = c323*h3 + c322*h2
  double c322,c323;
  c323 = dy1 * T32 * 2./(dx3+dx2);
  cout << "c323: " << c323 << endl;
  c322 = - c323;
  cout << "c322: " << c322 << endl;
  // Q52 = C525*h5 + C522*h2
  double c525,c522;
  c525 = dx2 * T52 * 2./(dy5+dy2);
  cout << "c525: " << c525 << endl;
  c522 = - c525;
  cout << "c522: " << c522 << endl;
  // C121*h1 + (C122+C322+C522)*h2 + C323*h3 + C525*h5 + QR + QP2 = 0
  // a21*h1 + a22*h2 + a23*h3 + a25*h5 + a20 = 0
  double a20,a21,a22,a23,a25;
  a21 = c121;
  a22 = c122+c322+c522;
  a23 = c323;
  a25 = c525;
  a20 = QR + QP2; // -b2
  double b20,b21,b23,b25;
  b21 = -a21/a22 * h1;
  cout << "b21: " << b21 << endl;
  b23 = -a23/a22;
  cout << "b23: " << b23 << endl;
  b25 = -a25/a22;
  cout << "b25: " << b25 << endl;
  b20 = -a20/a22;
  cout << "b20: " << b20 << endl;
  //===================================
  // Hausaufgabe an alle zum 18.05.2018
  // c522 = dx2*T52*2./dy2;
  // c525 = dx2*T52*2./dy1;
  //===================================
  // Hausaufgabe an alle zum 01.06.2018
  //-----------------------------------
  // Bilanz Knoten 3
  // 2.3: Q23+Q63+QR+QP3=0
  // Q23 = dy23 * T23 *(h2-h3)/(dx2/2+dx3/2)
  // Q23 = c232*h2 + c233*h3
  cout << "Knoten 3:" << endl;
  //...................................
  double c232,c233;
  double dy23 = dy1;
  double T23 = (dy2+dy3)/(dy2/T2+dy3/T3);
  c232 = dy23 * T23 / (dx2/2+dx3/2);
  cout << "c232: " << c232 << endl;
  c233 = - c232;
  cout << "c233: " << c233 << endl;
  // Q63 = dy63 * T63 *(h6-h3)/(dx6/2+dx3/2)
  // Q63 = c633*h3 + c636*h6
  double dx63 = dx3;
  double T63 = (dy2+dy6)/(dy2/T6+dy3/T3);
  double c633,c636;
  c633 = - dx63 * T63 / (dy6/2+dy3/2);
  cout << "c633: " << c633 << endl;
  c636 = - c633;
  //................
  double a30,a32,a33,a36;
  a30 = QR + QP3;
  a32 = c232;
  a33 = c233 + c633;
  a36 = c636;
  double b30,b32,b36;
  b30 = - a30 / a33;
  cout << "b30: " << b30 << endl;
  b32 = - a32 / a33;
  cout << "b32: " << b32 << endl;
  b36 = - a36 / a33;
  cout << "b36: " << b36 << endl;
  //-----------------------------------
  // Bilanz Knoten 5
  // 2.5: Q25+Q45+Q65+QR+QP5=0
  cout << "Knoten 5:" << endl;
  //...................................
  // Q25 = dx25 * T25 * (h2-h5)/(dy2/2+dy5/2)
  // Q25 = c252*h2 + c255*h5
  //...................................
  double c252,c255;
  double dx25 = dx2;
  double T25 = (dy2+dy5)/(dy2/T2+dy5/T5);
  c252 = dx25 * T25 / (dy2/2.+dy5/2.);
  cout << "c252: " << c252 << endl;
  c255 = - c252;
  cout << "c255: " << c255 << endl;
  //...................................
  // Q45 = dy45 * T45 * (h4-h5)/(dx4/2+dx5/2)
  // Q45 = c454*h4 + c455*h5
  //...................................
  double c454,c455;
  double dy45 = dy1;
  double T45 = (dx4+dx5)/(dx4/T4+dx5/T5);
  c454 = dy45 * T45 / (dx4/2.+dx5/2.);
  cout << "c454: " << c454 << endl;
  c455 = - c454;
  cout << "c455:" << c455 << endl;
  //...................................
  // Q65 = dy65 * T65 * (h6-h5)/(dx6/2+dx5/2)
  // Q65 = c656*h6 + c655*h5
  //...................................
  double c656,c655;
  double dy65 = dy1;
  double T65 = (dx6+dx5)/(dx6/T6+dx5/T5);
  c656 = dy65 * T65 / (dx6/2+dx5/2);
  cout << "c656: " << c656 << endl;
  c655 = - c656;
  cout << "c655:" << c655 << endl;
  //...................................
  double a50,a52,a54,a55,a56;
  a50 = QR + QP5;
  a52 = c252;
  a54 = c454;
  a55 = c255 + c455 + c655;
  a56 = c656;
  double b50,b52,b54,b56;
  b50 = - a50 / a55;
  cout << "b50: " << b50 << endl;
  b52 = - a52 / a55;
  cout << "b52: " << b52 << endl;
  b54 = - a54 / a55 * h4;
  cout << "b54: " << b54 << endl;
  b56 = - a56 / a55;
  cout << "b56: " << b56 << endl;
  //-----------------------------------
  // Bilanz Knoten 6
  // 2.6: Q36+Q56+QR+QP6=0
  // Q36 = dx36 * T36 *(h3-h6)/(dy3/2+dy6/2)
  // Q36 = c363*h3 + c366*h6
  cout << "Knoten 6:" << endl;
  //...................................
  double c363,c366;
  double dx36 = dx3;
  double T36 = (dy3+dy6)/(dy3/T3+dy6/T6);
  c363 = dx36 * T36 / (dy3/2+dy6/2);
  cout << "c363: " << c363 << endl;
  c366 = - c363;
  cout << "c366: " << c366 << endl;
  //...................................
  // Q56 = dy56 * T56 *(h5-h6)/(dx5/2+dx6/2)
  // Q56 = c565*h5 + c566*h6
  double dy56 = dy1;
  double T56 = (dx5+dx6)/(dx5/T5+dx6/T6);
  double c565,c566;
  c565 = dy56 * T56 / (dx5/2+dx6/2);
  cout << "c565: " << c565 << endl;
  c566 = - c565;
  //................
  double a60,a63,a65,a66;
  a60 = QR + QP6;
  a63 = c363;
  a65 = c565;
  a66 = c366 + c566;
  double b60,b63,b65;
  b60 = - a60 / a66;
  cout << "b60: " << b60 << endl;
  b63 = - a63 / a66;
  cout << "b63: " << b63 << endl;
  b65 = - a65 / a66;
  cout << "b65: " << b65 << endl;

  x[0]=15.;
  x[3]=15.;
  //...........................................................
  // Lösungsverfahren: Gauss-Seidel-Verfahren
    for(int k=0;k<solver_iterations;k++)
    {
      x0[1]=x[1];
      x0[2]=x[2];
      x0[4]=x[4];
      x0[5]=x[5];
      //x[1] = 0.2408 * x[2] + 0.3211 * x[4] + 4.4181;
      x[1] = b21 + b23 * x[2] + b25 * x[4] + b20;
      //x[2] = 0.4285 * x[1] + 0.5714 * x[5] + 0.0857;
      x[2] = b32 * x[1] + b36 * x[5] + b30;
      //x[4] = 0.7028 * x[1] + 0.1054 * x[5] + 2.0223 ;
      x[4] = b52 * x[1] + b54 + b56 * x[5] + b50;
      //x[5] = 0.8695 * x[2] + 0.1304 * x[4];
      x[5] = b63 * x[2] + b65 * x[4] + b60;
      // Ausgabe
      out_file << "Iteration step: " << k << endl;
      for(int i=0;i<n;i++)
      {
        out_file << "h" << i+1 << ":" << x[i] << endl;
      }
      //Fehlerberechung
      y = abs(x0[1]-x[1]) + abs(x0[2]-x[2]) + abs(x0[4]-x[4]) + abs(x0[5]-x[5]);
      y = y/4.;
      //cout << y << endl;
      if(y<eps) return 1;
    }
  //...........................................................
  // Output
  return 0;
}

