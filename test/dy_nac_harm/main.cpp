#include <iostream>
#include <qpbranch/dy_branch.hpp>
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
  int numI = 1;
  OpMat HeIJ(boost::extents[numI][numI]);
  OpMat XkIJP(boost::extents[numI][numI]);  
  VectorXd xs = VectorXd::LinSpaced(100, -5.0, 5.0);
  double k = 1.0;
  VectorXd ys = (k/2) * xs.array()*xs.array();  
  HeIJ[0][0] = new OperatorSpline(xs, ys);
  XkIJP[0][0] = nullptr;
  
  DyNac *dy = new DyNac(HeIJ, XkIJP, ns, "frozen", 2, 0.01);
  dy->m_ = 1.0;
  dy->q0_ = 1.0;
  dy->p0_ = 1.0;
  dy->SetUp();

  double dt = 0.3;
  int nt = 100;
  int n1t = 10;

  for(int it = 0; it < nt; it++) {
    con.write_f("t", it, it*dt);
    dy->DumpCon(it);
    for(int iit = 0; iit < n1t; iit++)
      dy->Update(dt/n1t);
  }
  
  delete dy;
  return 0;
}


