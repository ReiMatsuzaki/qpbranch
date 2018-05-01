/**
 *  cpp version of mangan4/con written by K.Y.
 */
#include <string>
#include <fstream>
#include <Eigen/Core>

namespace mangan4 {
  using std::string;
  using std::fstream;
  using Eigen::VectorXi;
  using Eigen::VectorXd;
  using Eigen::MatrixXi;
  using Eigen::MatrixXd;
  
  class Con {    
  private:    
    string root_;
    fstream file_;
    Con() : root_("_con") {}
    Con(const Con&);
    Con& operator=(const Con&);
    ~Con() {}
    string filepath(string key, int num);
  public:
    static Con& getInstance() {
      static Con singleton;
      return singleton;
    }
    void set_root(string root) { root_=root; }
    void open(string key, int num, char mode);
    void close();
    bool has(string key, int num) const;

    void write_l(string key, int num, bool val);
    void write_i(string key, int num, int val);
    void write_f(string key, int num, double val);
    void write_i1(string key, int num, const VectorXi& val);
    void write_i2(string key, int num, const MatrixXi& val);
    void write_f1(string key, int num, const VectorXd& val);
    void write_f2(string key, int num, const MatrixXd& val);
  };
}
