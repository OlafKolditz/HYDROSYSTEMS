#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <vector>
#include <fstream>
#include <QLineEdit>
#include <QRadioButton>
using namespace std;

class Dialog : public QDialog
{
    Q_OBJECT

public:
    Dialog(QWidget *parent = 0);
    ~Dialog();

private:
  int n;
  ofstream out_file;
  QPushButton* pushButtonRUN;
  double* A;
  double* b;
  double* x;
  int solver_method;
  int solver_iterations;
  QRadioButton* solver_method_;
  double eps;
private slots:
  void on_pushButtonRUN_clicked();
  void TestOutput(double*);
  void TestOutput(double*,double*);
  void GaussSeidel();
  void AssembleEQS();
  void CalculateFluxes();
};

#endif // DIALOG_H
