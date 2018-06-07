#include <iostream>
#include <qpbranch/dy_nac.hpp>
#include <qpbranch/con.hpp>

using namespace std;
using namespace qpbranch;
using mangan4::Con;

int main() {

  Con& con = Con::getInstance();

  // nuclear part
  int numA = 1;
  VectorXi ns(numA); ns << 0;

  // electron part
  int numI = 2;
  OpMat HeIJ(boost::extents[numI][numI]);
  OpMat XkIJP(boost::extents[numI][numI]);  
  VectorXd xs = VectorXd::LinSpaced(200, -10.0, 10.0);
  double A = 0.01;
  double B = 1.6;
  double C = 0.005;
  double D = 1.0;
  VectorXd y00s(xs.size());
  VectorXd y01s(xs.size());  
  VectorXd y11s(xs.size());
  VectorXd y10s(xs.size());
  for(int i = 0; i < xs.size(); i++) {
    double x = xs(i);
    if(x>0) 
      y00s(i) = A*(1-exp(-B*x));
    else
      y00s(i) = -A*(1-exp(+B*x));
    y01s(i) = C*exp(-D*x*x);
    y10s(i) = y01s(i);
    y11s(i) = -y00s(i);
  }
  HeIJ[0][0] = new OperatorSpline(xs, y00s);
  HeIJ[0][1] = new OperatorSpline(xs, y01s);
  HeIJ[1][0] = new OperatorSpline(xs, y10s);
  HeIJ[1][1] = new OperatorSpline(xs, y11s);
  XkIJP[0][0] = nullptr;
  XkIJP[0][1] = nullptr;
  XkIJP[1][1] = nullptr;
  XkIJP[1][0] = nullptr;

  // save
  con.write_f1("xs", 0, xs);

  // dynamics
  DyNac *dy = new DyNac(HeIJ, XkIJP, ns, "thawed", 2, 0.01);
  dy->type_intenuc_ = "RK4";
  dy->gamma0_ = 1.0;
  dy->q0_ = -7.0;
  dy->p0_ = +15.0;
  dy->m_  = 2000.0;
  dy->SetUp();

  // time grid
  double dt = 20.0;
  int nt = 100;
  int n1t = 1;

  // time loop start
  for(int it = 0; it < nt; it++) {
    con.write_f("t", it, it*dt);
    dy->DumpCon(it);
    for(int I = 0; I < 2; I++) {
      VectorXcd ys(xs.size());
      dy->At(I, xs, &ys);
      con.write_f1("psi_" + to_string(I) + "_re", it, ys.real());
      con.write_f1("psi_" + to_string(I) + "_im", it, ys.imag());
    }
    
    for(int iit = 0; iit < n1t; iit++) {
      dy->Update(dt/n1t);
    }
  }
  
  delete dy;
  return 0;
}


