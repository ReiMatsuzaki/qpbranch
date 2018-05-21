#include <complex>
#include <map>
#include <boost/multi_array.hpp>
#include <Eigen/Core>
#include "univariate.hpp"
#include "operator.hpp"

namespace qpbranch {

  using std::complex;
  using std::map;
  using boost::multi_array;
  using Eigen::MatrixXcd;
  using Eigen::VectorXi;

  typedef Univariate<complex<double> > Poly;

  // forward declaration
  class PolyPwgto;
  class OpBuffBasic;
  class OpBuffGausspot;

  // buffer for each operator
  class OpBuff {
  public:
    virtual int Maxn() const {return 0;}
    virtual void SetUp(PolyPwgto *) {}
    virtual void At( PolyPwgto*, const VectorXcd&, const VectorXd&, VectorXcd*) {}
    virtual void MatrixAsKet(OpBuff*, PolyPwgto*, MatrixXcd*) {}
    virtual void MatrixAsBra(OpBuffBasic*, PolyPwgto*, MatrixXcd*) {}
    virtual void MatrixAsBra(OpBuffGausspot*, PolyPwgto*, MatrixXcd*) {}
  };
  class OpBuffBasic : public OpBuff {
  protected:
    vector<Poly> poly_;
  public:
    OpBuffBasic(int num) : poly_(num) {}
    const Poly& poly(int A) { return poly_[A]; }
    int Maxn() const {return 0;}
    virtual void SetUp(PolyPwgto *basis) = 0;
    void At( PolyPwgto*, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);
    void MatrixAsKet(OpBuff *opbra, PolyPwgto *basis, MatrixXcd *res);
    void MatrixAsBra(OpBuffBasic *opket, PolyPwgto *basis, MatrixXcd *res);
    void MatrixAsBra(OpBuffGausspot *opket, PolyPwgto *basis, MatrixXcd *res);
  };
  // Buffer for id operator
  class OpBuffId : public OpBuffBasic {
  public:
    OpBuffId(int num) : OpBuffBasic(num) {}
    int Maxn() const { return 0; }
    void SetUp(PolyPwgto *basis);
  };
  // Buffer for Rn operator
  class OpBuffRn : public OpBuffBasic {
    int n_;
  public:    
    OpBuffRn(int num, int n) : OpBuffBasic(num), n_(n) {}
    int n() const { return n_; }
    int Maxn() const { return n_; }
    void SetUp(PolyPwgto *basis);
  };
  // Buffer for Pn operator
  class OpBuffPn : public OpBuffBasic {
    int n_;
  public:    
    OpBuffPn(int num, int n) : OpBuffBasic(num), n_(n) {}
    int n() const { return n_; }
    int Maxn() const { return n_; }
    void SetUp(PolyPwgto *basis);
  };
  // Buffer for Da operator
  class OpBuffDa : public OpBuffBasic {
    int id_;
  public:    
    OpBuffDa(int num, int id);    
    int id() const { return id_; }
    int Maxn() const;
    void SetUp(PolyPwgto *basis);
  };
  // Buffer for gauss potential
  class OpBuffGausspot : public OpBuff {
    OperatorGausspot *op_;
  public:
    OpBuffGausspot(OperatorGausspot *op) : op_(op) {}
    OperatorGausspot *op() const { return op_; }
    int Maxn() const { return 0; }
    void SetUp(PolyPwgto*) {}
    //    void MatrixAsKet(OpBuf *opbra, PolyPwgto *basis, MatrixXcd *res);
    //    void MatrixAsBra(OpBufBasic *opket, PolyPwgto *basis, MatrixXcd *res);
    //    void MatrixAsBra(OpBufGausspot *opket, PolyPwgto *basis, MatrixXcd *res);
  };  

  // Polynomial Plane Wave GTO basis set.
  // Each basis can be represented as
  //      GA = sum_i c_i (q-qA)^ni Exp[-gA(q-qA)^2 + i.pA(q-qA)]
  // for polynomial part Univariate class is used.
  class PolyPwgto {
  protected:    
    // data size
    int num_, nop_;
    // polynomial structure
    vector<Poly> polys_;
    vector<Operator*> ops_;
    // variable
    VectorXcd gs_;
    VectorXd Rs_, Ps_;
    // intermediate    
    bool is_setup_;
    Operator* op_id_;    
    map<Operator*, OpBuff*> buff_map_;
    VectorXd Ns_;
    VectorXi maxn_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
  public:
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
  public:
    PolyPwgto(const vector<Poly>& polys,
	      const vector<Operator*>& ops); 
    ~PolyPwgto();
    int num() const { return num_; }
    int nop() const { return nop_; }
    const vector<Operator*>& ops() const { return ops_; }
    const Poly& poly(int A) const { return polys_[A]; }
    const VectorXcd& gs() const { return gs_; }
    const VectorXd& Rs() const { return Rs_; }
    const VectorXd& Ps() const { return Ps_; }
    const VectorXd& Ns() const { return Ns_; }
    VectorXcd& ref_gs() { return gs_; }
    void SetUp();
    int size() const { return num_; }
    void Matrix(Operator *ibra, Operator *iket, MatrixXcd *res);
    void At(Operator *iop, const VectorXcd& cs, const VectorXd& xs, VectorXcd *res);    
    void DumpCon(int it);
    inline complex<double> get_d(int A, int B,int na,int nb,int Nk) {
      return (*(*d_)[A][B])[na][nb][Nk];
    }
  };

  
}
