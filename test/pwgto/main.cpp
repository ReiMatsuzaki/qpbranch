#include <iostream>

#include <gtest/gtest.h>

#include <qpbranch/gtestplus.hpp>
#include <qpbranch/mathplus.hpp>
#include <qpbranch/eigenplus.hpp>
#include <qpbranch/pwgto.hpp>
#include <qpbranch/pwgto1c.hpp>
#include <qpbranch/pwgto_buf.hpp>
#include <qpbranch/ppwgto.hpp>

using namespace std;
using namespace qpbranch;
using namespace boost;

TEST(utest_pwgto, test_coef_d) {
  complex<double> gAB(1, 0.2);
  complex<double> wAB(1.2, 0);
  double RA(0.3), RB(0.3);
  const int maxnA = 3;
  const int maxnB = 2;
  multi_array<complex<double>, 3> d(extents[maxnA+1][maxnB+1][maxnA+maxnB+1]);;
  HermiteCoefD(gAB, wAB, RA, RB, maxnA, maxnB, &d);
  
  for(int nA=0; nA<=maxnA; nA++) {
    for(int nB=0; nB<=maxnB; nB++) {
      for(int Nk=0; Nk<=maxnA+maxnB; Nk++) {
	auto ref = HermiteCoefD0(gAB, wAB, RA, RB, nA, nB, Nk);
	ASSERT_DOUBLE_EQ(real(ref), real(d[nA][nB][Nk])) << "nA:" << nA << endl;
	ASSERT_DOUBLE_EQ(imag(ref), imag(d[nA][nB][Nk])) << "nA:" << nA << endl;	
      }
    }
  }

}
TEST(utest_pwgto, overlap) {

  int num = 4;
  VectorXi ns(num); ns << 0, 0, 2, 2;
  //  OperatorId opopid;  
  //  auto *op_id = &opopid;
  auto *op_id = new OperatorId();
  
  vector<Operator*> ops = {op_id};
  auto *basis = new Pwgto(ns, ops);
  basis->ref_gs() << 1.0,  complex<double>(0.9,-0.8), 1.0, complex<double>(0.2, 0.1);
  basis->ref_Rs() << 0.0, 0.0, 0.2, 0.2;
  basis->ref_Ps() << 0.0, 0.0, 0.5, 0.5;
  basis->SetUp();

  MatrixXcd S(num, num);
  basis->Matrix(op_id, op_id, &S);

  for(int A = 0; A<num; A++) {
    ASSERT_DOUBLE_EQ(1.0, S(A,A).real()) << "A:" << A << endl;
    ASSERT_DOUBLE_EQ(0.0, S(A,A).imag()) << "A:" << A << endl;
  }

  for(int A = 0; A<num; A++) {
    for(int B = 0; B<A; B++) {
      ASSERT_DOUBLE_EQ(real(S(A,B)), +real(S(B,A))) << "A:" << A << endl << "B:" << B << endl;
      ASSERT_DOUBLE_EQ(imag(S(A,B)), -imag(S(B,A))) << "A:" << A << endl << "B:" << B << endl;
    }
  }

  int A = 2;
  int B = 3;
  auto calc = S(A, B)/(basis->Ns()(A)*basis->Ns()(B));
  VectorXcd gg(9+1);
  IntGto2N(9, conj(basis->gs()(A))+basis->gs()(B), &gg);
  auto ref = gg((basis->ns()(A) + basis->ns()(B))/2);

  ASSERT_DOUBLE_EQ(real(ref), real(calc)) << "A:" << A << endl << "B:" << B << endl;
  ASSERT_DOUBLE_EQ(imag(ref), imag(calc)) << "A:" << A << endl << "B:" << B << endl;
  
}
TEST(utest_pwgto, multipole) {

  int num = 5;
  auto id = new OperatorId();
  auto R1 = new OperatorRn(1);
  auto R2 = new OperatorRn(2);
  vector<Operator*> ops = {id, R1, R2};
  VectorXi ns(num); ns << 0, 0, 2, 1, 0;
  auto *basis = new Pwgto(ns, ops);
  basis->ref_gs() << 1.1, 1.1, 1.2, 0.4, 0.8;
  basis->ref_Rs() << 0.0, 0.1, 0.1, 0.1, 0.3;  
  basis->ref_Ps() << 0.0, 0.0, 0.3, 0.0, 0.1;
  basis->SetUp();

  MatrixXcd S(num, num);
  basis->Matrix(id, id, &S);
  
  MatrixXcd M01(num, num);
  basis->Matrix(id, R1, &M01);

  MatrixXcd M10(num, num);
  basis->Matrix(R1, id, &M10);

  MatrixXcd M02(num, num);
  basis->Matrix(id, R2, &M02);

  MatrixXcd M11(num, num);
  basis->Matrix(R1, R1, &M11);

  MatrixXcd M20(num, num);
  basis->Matrix(R2, id, &M20);
  
  for(int A = 0; A < num; A++) {
    auto g = conj(basis->gs()(A)) + basis->gs()(A);
    int n = basis->ns()(A);
    VectorXcd gg(n+1+1);
    IntGto2N(n+1, g, &gg);
    auto ref = gg(n+1) * pow(basis->Ns()(A), 2);
    ASSERT_DOUBLE_EQ(real(ref), real(M20(A,A)));
    ASSERT_DOUBLE_EQ(real(ref), real(M02(A,A)));
  }

  int A = 1;
  int B = 2;
  ASSERT_DOUBLE_EQ(real(M01(A,B)), real(M10(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  ASSERT_DOUBLE_EQ(imag(M01(A,B)), imag(M10(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  
  ASSERT_DOUBLE_EQ(real(M02(A,B)), real(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  ASSERT_DOUBLE_EQ(imag(M02(A,B)), imag(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  
  ASSERT_DOUBLE_EQ(real(M20(A,B)), real(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;
  ASSERT_DOUBLE_EQ(imag(M20(A,B)), imag(M11(A,B))) << "A:"<<A<<endl<<"B:"<<B<<endl;

  delete basis;
  
}
TEST(utest_pwgto, pn) {
  
  int num = 5;
  auto id = new OperatorId();
  auto P1 = new OperatorPn(1);
  auto P2 = new OperatorPn(2);
  vector<Operator*> ops = {id, P1, P2};
  VectorXi ns(num); ns << 0, 0, 2, 1, 0;
  auto *basis = new Pwgto(ns, ops);
  basis->ref_gs() << 1.1, 1.1, 1.2, 0.4, 0.8;
  basis->ref_Rs() << 0.0, 0.1, 0.1, 0.1, 0.3;  
  basis->ref_Ps() << 0.0, 0.0, 0.3, 0.0, 0.1;
  basis->SetUp();

  MatrixXcd M01(num, num);
  basis->Matrix(id, P1, &M01);

  MatrixXcd M10(num, num);
  basis->Matrix(P1, id, &M10);

  MatrixXcd M02(num, num);
  basis->Matrix(id, P2, &M02);

  MatrixXcd M11(num, num);
  basis->Matrix(P1, P1, &M11);

  MatrixXcd M20(num, num);
  basis->Matrix(P2, id, &M20);
  
  for(int A = 0; A < num; A++) {
    for(int B = 0; B < num; B++) {
      double tol = pow(10.0, -10);
      ASSERT_NEAR(real(M01(A,B)), real(M10(A,B)), tol) <<"A:"<<A<<endl<<"B:"<<B;
      ASSERT_NEAR(real(M11(A,B)), real(M02(A,B)), tol) <<"A:"<<A<<endl<<"B:"<<B;
      ASSERT_NEAR(real(M11(A,B)), real(M20(A,B)), tol) <<"A:"<<A<<endl<<"B:"<<B;
    }
  }

  delete basis;
}
TEST(utest_pwgto, da) {

  int num = 4;
  
  complex<double> b(1.3);
  complex<double> v0(1.2);
  complex<double> q0(0.0);    

  double dx = 0.01;
  double R0 = 0.1;
  double R1 = 0.2;
  double P0 = 0.3;
  double P1 = 0.2;
  complex<double> g0 = complex<double>(1.0, 0.1);
  complex<double> g1 = complex<double>(1.3, 0.15);
  
  for(int n = 0; n < 2; n++) {
    VectorXi ns = VectorXi::Constant(num, n);
    vector<int> ops_id = {kIdDR, kIdDP, kIdDgr, kIdDgi};  
    for(auto it = ops_id.begin(); it!=ops_id.end(); ++it) {
      
      int kId = *it;
      auto op_id = new OperatorId();
      auto op_da = new OperatorDa(kId);
      auto op_pot = new OperatorGausspot(v0, b, q0);
      vector<Operator*> ops = {op_id, op_da, op_pot};
      
      auto *basis = new Pwgto(ns, ops);
      basis->ref_Rs() << R0, R0, R0, R1;
      basis->ref_Ps() << P0, P0, P0, P1;
      basis->ref_gs() << g0, g0, g0, g1;
      if(kId==kIdDR) {
	basis->ref_Rs()(1) += dx;
	basis->ref_Rs()(2) -= dx;
      } else if(kId==kIdDP) {
	basis->ref_Ps()(1) += dx;
	basis->ref_Ps()(2) -= dx;
      } else if(kId==kIdDgr) {
	basis->ref_gs()(1) += complex<double>(dx, 0.0);
	basis->ref_gs()(2) -= complex<double>(dx, 0.0);
      } else if(kId==kIdDgi) {
	basis->ref_gs()(1) += complex<double>(0.0, dx);
	basis->ref_gs()(2) -= complex<double>(0.0, dx);
      }
      basis->SetUp();
      
      MatrixXcd M1(num,num), M2(num,num);
      basis->Matrix(op_id,  op_da, &M1);
      basis->Matrix(op_id,  op_id, &M2);
      EXPECT_NEAR(M1(3, 0).real(), ((M2(3,1)-M2(3,2))/(2*dx)).real(), dx*dx) <<
	"kId:" << kId << endl <<
	"n  :" << n   << endl;
      
      basis->Matrix(op_da, op_pot, &M1);
      basis->Matrix(op_id, op_pot, &M2);
      EXPECT_NEAR(M1(0, 3).real(), ((M2(1,3)-M2(2,3))/(2*dx)).real(), dx*dx) <<
	"kId:" << kId << endl <<
	"n  :" << n   << endl;
    }
  }
}
TEST(utest_pwgto, test_harmonic) {

  int num = 4;
  double x0 = 0.3d;
  double k = 0.39;
  double m = 2.0;
  double w = sqrt(k/m);

  VectorXi ns(num); ns << 0, 1, 2, 3;
  auto opid = new OperatorId();
  auto opR2 = new OperatorRn(2);
  auto opP2 = new OperatorPn(2);  
  vector<Operator*> ops = {opid, opR2, opP2};
  auto *basis = new Pwgto(ns, ops);

  complex<double> g = m*w/2;
  basis->ref_gs() = VectorXcd::Ones(num)*g;
  basis->ref_Rs() = VectorXd::Ones(num )*x0;
  basis->ref_Ps() = VectorXd::Zero(num);
  basis->SetUp();

  MatrixXcd P2(num,num), R2(num,num), H(num,num), S(num,num);
  basis->Matrix(opid, opid, &S);
  basis->Matrix(opid, opR2, &R2);
  basis->Matrix(opid, opP2, &P2);
  H = P2/(2*m) + k/2*R2;

  VectorXd eigs;
  MatrixXcd U;
  Zhegv(H, S, &U, &eigs);

  int nx(3);
  VectorXd xs(nx);
  xs << 0.2, 0.25, 0.3;

  for(int n = 0; n < num-1; n++) {
    ASSERT_DOUBLE_EQ(w*(n+0.5), eigs(n)) << "n: " << n;
    VectorXcd ys(nx);
    basis->At(opid, U.col(n), xs, &ys);
  }
  
  delete basis;  
  
}
TEST(utest_pwgto, test_gausspot) {

  int num = 1;
  complex<double> b(1.3);
  complex<double> v0(1.2);
  complex<double> q0(0.0);
  
  VectorXi ns(num); ns << 0;
  auto opid = new OperatorId();
  auto opdR = new OperatorDa(kIdDR);
  auto opV  = new OperatorGausspot(v0, b, q0);
  vector<Operator*> ops = {opid, opdR, opV};
  auto *basis = new Pwgto(ns, ops);
  basis->ref_gs() << 1000.0;
  basis->ref_Rs() << 0.3;
  basis->SetUp();
  
  MatrixXcd V(num,num);
  basis->Matrix(opid, opV, &V);
  double R0 = basis->Rs()(0);

  double dR = 0.001;
  VectorXd Rs(3); Rs << R0, R0-dR, R0+dR;
  VectorXcd refs(3);
  opV->At(Rs, &refs);
  ASSERT_NEAR(real(refs(0)), real(V(0,0)), 3.0*pow(10.0, -3));

  basis->Matrix(opdR, opV, &V);
  complex<double> ref = (refs[2]-refs[1])/(2.0*dR);
  ASSERT_NEAR(real(ref), 2.0*real(V(0,0)), pow(10.0, -3));

  delete basis;
  
}
TEST(utest_pwgto, one_center) {
  
  complex<double> b(1.3);
  complex<double> v0(1.2);
  complex<double> q0(0.0);
  auto opid = new OperatorId();
  auto opR1 = new OperatorRn(1);
  auto opR2 = new OperatorRn(2);
  auto opP1 = new OperatorPn(1);
  auto opP2 = new OperatorPn(2);
  auto opDR = new OperatorDa(kIdDR);
  auto opDP = new OperatorDa(kIdDP);
  auto opDgr= new OperatorDa(kIdDgr);
  auto opDgi= new OperatorDa(kIdDgi);  
  auto opV  = new OperatorGausspot(v0, b, q0);
  vector<Operator*> ops = {opid, opR1, opR2, opP1, opP2, opDR, opDP, opDgr, opDgi, opV};

  double R0 = 0.1;
  double P0 = 2.0;
  complex<double> g0(1.0,0.2);
  int num = 3;
  VectorXi ns(num); ns << 0, 1, 2;
  
  auto *pwgto   = new Pwgto(ns, ops);
  pwgto->ref_Rs() = VectorXd::Constant( num, R0);
  pwgto->ref_Ps() = VectorXd::Constant( num, P0);
  pwgto->ref_gs() = VectorXcd::Constant(num, g0);
  pwgto->SetUp();
  
  auto *pwgto1c = new Pwgto1c(  ns, R0, P0, g0, ops);
  pwgto1c->SetUp();

  double tol = pow(10.0, -12.0);

  MatrixXcd X(num,num);
  MatrixXcd Y(num,num);
  for(auto op1 : ops) {
    for(auto op2 : ops) {
      if(op1==opV)
	continue;
      pwgto->Matrix(op1, op2, &X);
      pwgto1c->Matrix(op1, op2, &Y);
      for(int A = 0; A < num; A++) {
	for(int B = 0; B < num; B++) {
	  ASSERT_NEAR(X(A,B).real(), Y(A,B).real(), tol) <<
	    "op1:" << op1->str() << endl <<
	    "op2:" << op2->str() << endl <<
	    " A :" << A << endl <<
	    " B :" << B << endl;
	  ASSERT_NEAR(X(A,B).imag(), Y(A,B).imag(), tol) <<
	    "op1:" << op1->str() << endl <<
	    "op2:" << op2->str() << endl <<
	    " A :" << A << endl <<
	    " B :" << B << endl;
	}
      }
    }
  }

  int nx = 2;
  VectorXd  xs(nx); xs << 0.1, 0.2;
  VectorXcd y1s(nx);
  VectorXcd y2s(nx);
  VectorXcd cs(num); cs << 0.5, complex<double>(0.4, 0.1), 0.2;
  for(auto op : ops) {
    if(op==opV)
      continue;
    pwgto->At(  op, cs, xs, &y1s);
    pwgto1c->At(op, cs, xs, &y2s);
    for(int i = 0; i < nx; i++) {
      ASSERT_DOUBLE_EQ(y1s(i).real(), y2s(i).real()) <<
	"op:" << op->str() << endl <<
	" i :"   << i << endl;
      ASSERT_DOUBLE_EQ(y1s(i).imag(), y2s(i).imag()) <<
	"op:" << op->str() << endl <<
	" i :"   << i << endl;
    }
  }
  
}
TEST(utest_pwgto, test_poly) {

  /*
  
  complex<double> b(1.3);
  complex<double> v0(1.2);
  complex<double> q0(0.0);
  auto opid = new OperatorId();
  auto opR1 = new OperatorRn(1);
  auto opR2 = new OperatorRn(2);
  auto opP1 = new OperatorPn(1);
  auto opP2 = new OperatorPn(2);
  auto opDR = new OperatorDa(kIdDR);
  auto opDP = new OperatorDa(kIdDP);
  auto opV  = new OperatorGausspot(v0, b, q0);
  vector<Operator*> ops = {opid, opR1, opR2, opP1, opP2, opDR, opDP, opV};

  double R0 = 0.1;
  double P0 = 2.0;
  complex<double> g0(1.0,0.2);
  int num = 3;
  VectorXi ns(num); ns << 0, 1, 2;
    
  auto *basis1 = new Pwgto1c(  ns, R0, P0, g0, ops);
  basis1->SetUp();

  vector<Poly> polys;
  polys.push_back(Poly({{1.0, 0}}));
  polys.push_back(Poly({{1.0, 1}}));
  polys.push_back(Poly({{1.0, 2}}));
  auto *basis2 = new PolyPwgto(polys, ops);
  basis2->SetUp();
  
  double tol = pow(10.0, -12.0);

  MatrixXcd M1(num,num);
  MatrixXcd M2(num,num);
  for(auto op1 : ops) {
    for(auto op2 : ops) {
      if(op1==opV)
	continue;
      basis1->Matrix(op1, op2, &M1);
      basis2->Matrix(op1, op2, &M2);
      for(int A = 0; A < num; A++) {
	for(int B = 0; B < num; B++) {
	  ASSERT_NEAR(M1(A,B).real(), M2(A,B).real(), tol) <<
	    "op1:" << op1->str() << endl <<
	    "op2:" << op2->str() << endl <<
	    " A :" << A << endl <<
	    " B :" << B << endl;
	  ASSERT_NEAR(M1(A,B).imag(), M2(A,B).imag(), tol) <<
	    "op1:" << op1->str() << endl <<
	    "op2:" << op2->str() << endl <<
	    " A :" << A << endl <<
	    " B :" << B << endl;
	}
      }
    }
  }

  */
  
}
TEST(utest_pwgto, test_spline) {

  auto k = 1.0;
  VectorXd xs = VectorXd::LinSpaced(100, -5.0, 5.0);
  VectorXd ys = k/2*(xs.array()*xs.array());
  auto op_v = new OperatorSpline(xs, ys);

  auto id = new OperatorId();
  auto R2 = new OperatorRn(2);
  vector<Operator*> ops = {id, op_v, R2};

  int num(2);
  VectorXi ns(num); ns << 0, 2;
  
  auto *basis = new Pwgto(ns, ops, 0, 0.01);  
  basis->ref_gs() << 1.0, 1.2;
  basis->ref_Rs() << 0.0, 0.0;
  basis->ref_Ps() << 0.0, 0.0;
  basis->SetUp();

  MatrixXcd M1(num, num);
  basis->Matrix(id, op_v, &M1);
  EXPECT_NEAR(abs(M1(0,0)), 0.0, pow(10.0, -12.0));

  auto *basis2 = new Pwgto(ns, ops, 2, 0.01);  
  basis2->ref_gs() << 1.0, 1.2;
  basis2->ref_Rs() << 0.0, 0.0;
  basis2->ref_Ps() << 0.0, 0.0;
  basis2->SetUp();

  MatrixXcd M2(num, num);
  basis2->Matrix(id, op_v, &M1);
  basis2->Matrix(id, R2,   &M2);
  M2 *= (k/2);

  EXPECT_MATRIXXCD_NEAR(M1, M2, pow(10.0, -10.0));
}
TEST(utest_pwgto, test_spline_dx) {

  auto k = 1.0;
  VectorXd xs = VectorXd::LinSpaced(100, -5.0, 5.0);
  VectorXd ys = k/2*(xs.array()*xs.array());
  auto op_v    = new OperatorSpline(xs, ys);
  auto op_v_p1 = new OperatorSplineP1(xs, ys);

  auto id = new OperatorId();
  auto p1 = new OperatorPn(1);
  vector<Operator*> ops = {id, p1, op_v, op_v_p1};

  int num(2);
  VectorXi ns(num); ns << 0,1;

  int norder = 0;
  auto *basis = new Pwgto(ns, ops, norder, 0.01);
  basis->ref_gs() << 1.0, 1.2;
  basis->ref_Rs() << 0.3, 0.0;
  basis->ref_Ps() << 0.3, 0.0;
  basis->SetUp();
  
  /*
  cout << endl << endl;
  cout << "id" << endl;
  basis->buffer_map_[id]->Dump();
  cout << endl << "op_v" << endl;
  basis->buffer_map_[op_v]->Dump();
  cout << endl << "op_p1" << endl;
  basis->buffer_map_[p1]->Dump();
  cout << endl << "op_v_p1" << endl;
  basis->buffer_map_[op_v_p1]->Dump();
  */
  
  MatrixXcd M1(num, num);
  basis->Matrix(id, op_v_p1, &M1);
  MatrixXcd M2(num, num);
  basis->Matrix(p1, op_v, &M2);
  EXPECT_MATRIXXCD_EQ(M1, M2);
  
  basis->Matrix(id, p1, &M1);
  basis->Matrix(p1, id, &M2);
  EXPECT_MATRIXXCD_EQ(M1, M2);
  
  delete basis;
    
}
