#include "dialog.h"

#include <QVBoxLayout>
#include <QPushButton>
#include <QRadioButton>
#include <QGroupBox>
#include <QMessageBox>
#include <QLabel>

#include <string>
using namespace std;
#include <cmath>
#include <iostream>

extern void Gauss(double*, double*, double*, int);

Dialog::Dialog(QWidget *parent) : QDialog(parent)
{
  //elements
  QLabel* labelHeader = new QLabel(tr("Picture source: Sebastian Bauer CAU Kiel"));
  QLabel *label_ogs = new QLabel();
  label_ogs->setAlignment(Qt::AlignCenter);
  label_ogs->setPixmap(QPixmap("ogs_teaching_150.png"));
  QLabel *label_exercise = new QLabel();
  label_exercise->setAlignment(Qt::AlignCenter);
  label_exercise->setPixmap(QPixmap("gw1_300.png"));
  pushButtonRUN = new QPushButton(tr("Run simulation"));
  solver_method_ = new QRadioButton("Solver method",parent);
  QGroupBox *groupBox = new QGroupBox(tr("Exclusive Radio Buttons"));
  QRadioButton *radio1 = new QRadioButton(tr("&Gauss Solver"));
  QRadioButton *radio2 = new QRadioButton(tr("Gauss-&Seidel Solver"));
  QRadioButton *radio3 = new QRadioButton(tr("&Alternative Solver"));
  QRadioButton *radio4 = new QRadioButton(tr("&Neues Verfahren"));
  radio2->setChecked(true);
  //connect
  connect(pushButtonRUN,SIGNAL(clicked()),this,SLOT(on_pushButtonRUN_clicked()));
  connect(radio1,SIGNAL(clicked(bool)),this,SLOT(clickkedstate(bool)));
  //layout
  QVBoxLayout *mainLayout = new QVBoxLayout;
  QVBoxLayout *upperLayout = new QVBoxLayout;
  QHBoxLayout *lowerLayout = new QHBoxLayout; //H->V
  QVBoxLayout *lowerleftLayout = new QVBoxLayout;
  upperLayout->addWidget(labelHeader);
  upperLayout->addWidget(label_exercise);
  upperLayout->addWidget(pushButtonRUN);
  lowerLayout->addLayout(lowerleftLayout);
  lowerleftLayout->addWidget(radio1);
  lowerleftLayout->addWidget(radio2);
  lowerleftLayout->addWidget(radio3);
  lowerleftLayout->addWidget(radio4);
  lowerleftLayout->addStretch(1);
  groupBox->setLayout(lowerleftLayout);
  lowerLayout->addWidget(label_ogs);
  mainLayout->addLayout(upperLayout);
  mainLayout->addLayout(lowerLayout);
  setLayout(mainLayout);
  //initializations I
  solver_method = 1;
  solver_iterations = 100;
  n=6;
  int n2=n*n;
  out_file.open("out.txt");
  //out_file.setf(ios::scientific);
  out_file.precision(5);
  eps = 1e-3;
  //data structures
  A = new double[n2]; //matrix
  x = new double[n]; //solution vector
  b = new double[n]; //RHS vector
  //initializations II
  for(int i=0;i<n;i++)
    x[i]=10.; //?
  x[0]=10.;
  x[3]=10.;
}

Dialog::~Dialog()
{
    delete [] A;
}

void Dialog::on_pushButtonRUN_clicked()
{
  QMessageBox msgBox;
  switch(solver_method)
  {
    case 0: //Gauss
      AssembleEQS();   //assemble equation system
      TestOutput(A,b);
      Gauss(A,b,x,n);  //solve EQS via Gauss
      break;
    case 1: //Gauss-Seidel
      //CalculateFluxes();
      GaussSeidel();
      msgBox.setText("Gauss-Seidel method finished, \n results in out.txt");
      break;
    case 3: // neues Verfahren
      msgBox.setText("Neues Verfahren vorbereitet");
     break;
  }
  pushButtonRUN->setStyleSheet("background-color: green");
  msgBox.exec();
}

void Dialog::GaussSeidel()
{
  double x0[6];
  double y;
  for(int k=0;k<solver_iterations;k++)
  {
    x0[1]=x[1];
    x0[2]=x[2];
    x0[4]=x[4];
    x0[5]=x[5];
    x[1] = 0.2408 * x[2] + 0.3211 * x[4] + 4.4181;
    x[2] = 0.4285 * x[1] + 0.5714 * x[5] + 0.0857;
    x[4] = 0.7028 * x[1] + 0.1054 * x[5] + 2.0223 ;
    x[5] = 0.8695 * x[2] + 0.1304 * x[4];
    out_file << "Iteration step: " << k << endl;
    TestOutput(x);
    //Fehlerberechung
    y = abs(x0[1]-x[1]) + abs(x0[2]-x[2]) + abs(x0[4]-x[4]) + abs(x0[5]-x[5]);
    y = y/4.;
    cout << y << endl;
    if(y<eps) return;
  }
}

void Dialog::TestOutput(double*x)
{
  for(int i=0;i<n;i++)
  {
    out_file << "h" << QString::number(i).toStdString() << ":" << x[i] << endl;
  }
}

void Dialog::TestOutput(double*a,double*b)
{

  for(int i=0;i<n;i++)
  {
    for(int j=0;j<n;j++)
    {
      out_file << A[n*i+j] << "\t";
    }
    out_file << " : " << b[i] << endl;
  }
}

void Dialog::AssembleEQS()
{
}

void Dialog::CalculateFluxes()
{
  //defintions
  float dx1,dx2,dx3,dy1,dy2,T1,T2,h1,h2,Q12;
  //values
  dx1 = 100.;
  dx2 = 1000.;
  dx3 = 1000.;
  dy1 = 500.;
  dy2 = 500.;
  T1 = 0.01;
  T2 = 0.05;
  h1 = 10.;
  h2 = 11.;
  Q12 = dy1*((dx1+dx2)/(dx1/T1+dx2/T2))*((h1-h2)/(dx1/2.+dx2/2.));
  cout << Q12 << endl;
}
