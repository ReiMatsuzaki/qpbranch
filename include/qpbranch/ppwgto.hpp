namespace qpbranch {

  // Polynomial Plane Wave GTO basis set.
  // Each basis can be represented as
  //      GA = sum_i c_i (q-qA)^ni Exp[-gA(q-qA)^2 + i.pA(q-qA)]
  class PolyPwgto {
  protected:
    // data size
    int num_, nop_;
    // polynomial structure
    typedef vector<vector<pair<int,complex<double> > > > Polys;
    Polys polys_;
    // variable
    VectorXcd gs_;
    // intermediate
    bool is_setup_;
    map<Operator*, OpBuf*> buf_map_;
    VectorXd Ns_;
    VectorXi maxn_;
    multi_array<multi_array<complex<double>,3>*,2> *d_;
  public:
    MatrixXcd gAB_, eAB_, hAB_,  RAB_;
  public:
    PolyPwgto(const Polys& polys, const vector<Operator*>& ops); 
    ~PolyPwgto();
    int num() const { return num_; }
    int nop() const { return nop_; }
    const vector<Operator*>& ops() const { return ops_; }
    OpBuf *poly() const { return poly; }
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
