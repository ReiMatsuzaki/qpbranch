#include <iostream>
#include <gtest/gtest.h>
#include <qpbranch/dy_mono.hpp>
#include <qpbranch/con.hpp>

using namespace std;
using namespace qpbranch;
using namespace mangan4;

TEST(TestDyMono, FreeParticle) {

  Con& con = Con::getInstance();
  con.set_root("_con_fp");

  // basic info
  int num = 1;
  VectorXi ns(num); ns << 0;

  // system (Free particle)
  auto m = 2.0;
  VectorXd xs = VectorXd::LinSpaced(400, -10.0, 10.0);
  VectorXd ys = 0.0*xs;
  auto vop = new OperatorSpline(xs, ys);

  // initial value
  auto q0 = -2.0;
  auto p0 = +1.0;
  auto g0 = complex<double>(1.2, 0.0);
  
  // dynamics
  auto dy = new DySetPoly(vop, ns, "thawed", 2, 0.001);
  dy->type_intenuc_ = "RK4";
  dy->gr0_ = real(g0);
  dy->gi0_ = imag(g0);
  dy->q0_  = q0;
  dy->p0_  = p0;
  dy->m_   = m;
  dy->SetUp();

  // time
  auto dt = 0.01;
  auto nt = 100;

  // exact
  double p, q;
  complex<double> g;

  // propagation
  for(int it = 0; it < nt; it++) {

    // time
    con.write_f("t", it, it*dt);
    
    // analytic solution. see Tannor's book p.25
    auto t = dt*it;
    auto i = complex<double>(0,1);
    g = g0 / (1.0+(2.0*i*g0*t/m));
    q = q0+p0*t/m;
    p = p0;
    con.write_f("exact_q0", it, q);
    con.write_f("exact_p0", it, p);
    con.write_f("exact_gr0", it, real(g));
    con.write_f("exact_gi0", it, imag(g));
    
    // dump calculation
    dy->DumpCon(it, "");    

    // update
    if(it!=nt-1)
      dy->Update(dt);
  }

  // compare
  auto tol = 0.00002;
  EXPECT_NEAR(real(g), dy->gr0_, tol);
  EXPECT_NEAR(imag(g), dy->gi0_, tol);
  EXPECT_NEAR(q,       dy->q0_,  tol);
  EXPECT_NEAR(p,       dy->p0_,  tol);
  
  delete dy;
}
TEST(TestDyMono, Harmonic) {

  Con& con = Con::getInstance();
  con.set_root("_con_harm");
  auto i = complex<double>(0,1);
  
  // basic info
  int num = 1;
  VectorXi ns(num); ns << 0;

  // system (Harmonics)
  auto m = 2.0;
  auto w = 0.8;  
  auto a = m*w/2;
  VectorXd xs = VectorXd::LinSpaced(400, -10.0, 10.0);
  VectorXd ys = (m*w*w/2)*xs.array()*xs.array();
  auto vop = new OperatorSpline(xs, ys);

  // initial value
  auto q0 = -2.0;
  auto p0 = +1.0;
  auto g0 = complex<double>(1.2, 0.0);
  
  // dynamics
  auto dy = new DySetPoly(vop, ns, "thawed", 2, 0.001);
  dy->type_intenuc_ = "RK4";
  dy->gr0_ = real(g0);
  dy->gi0_ = imag(g0);
  dy->q0_  = q0;
  dy->p0_  = p0;
  dy->m_   = m;
  dy->SetUp();

  // time
  auto dt = 0.1;
  auto nt = 100;
  
  // exact
  double q, p;
  complex<double> g;

  // propagation
  for(int it = 0; it < nt; it++) {

    // time
    con.write_f("t", it, it*dt);

    // analytic solutions. see Tannor's book p.29
    auto t = it*dt;
    auto cwt = cos(w*t);
    auto swt = sin(w*t);     
    g = a * (g0*cwt+i*a*swt) / (i*g0*swt+a*cwt);
    q = q0*cwt + p0/(m*w)*swt;
    p = p0*cwt - (m*w*q0)*swt;

    // dump calculation results
    dy->DumpCon(it, "");

    // dump exact
    con.write_f("exact_q0", it, q);
    con.write_f("exact_p0", it, p);
    con.write_f("exact_gr0", it, real(g));
    con.write_f("exact_gi0", it, imag(g));

    // update
    if(it!=nt-1)
      dy->Update(dt);
  }

  // compare
  auto tol = 0.00002;
  EXPECT_NEAR(real(g), dy->gr0_, tol);
  EXPECT_NEAR(imag(g), dy->gi0_, tol);
  EXPECT_NEAR(q,       dy->q0_,  tol);
  EXPECT_NEAR(p,       dy->p0_,  tol);
  
  delete dy;  
}

