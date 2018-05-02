#include <sys/stat.h>
#include <iostream>

#include <boost/format.hpp>
#include <qpbranch/con.hpp>

namespace mangan4 {
  using namespace std;
  using namespace Eigen;

  template<class Scalar>
  void write_vector(fstream& ofs, const Matrix<Scalar,Dynamic,1>& val) {
    ofs << "[";
    for(int i = 0; i < val.size(); i++) {
      if(i == 0) 
	ofs << " " << endl;
      else 
	ofs << ", ";
      ofs << val[i];
    }    
    ofs << endl << "]" << endl;
  }
  template<class Scalar>
  void write_matrix(fstream& ofs, const Matrix<Scalar,Dynamic,Dynamic>& val) {
    
    ofs << "[";
    for(int i = 0; i < val.cols(); i++) {      
      if(i == 0)
	ofs << " " << endl;
      else 
	ofs << ", ";

      /*
      Matrix<Scalar,Dynamic,1>& vec = val.col(i);
      write_vector(ofs, vec);
      */

      ofs << "[";
      for(int j = 0; j < val.rows(); j++) {
	
	if(j == 0)
	  ofs << " " << endl;
	else
	  ofs << ", ";
	ofs << val(j,i);
      }
      ofs << endl << "]" << endl;

      
    }
    ofs << "]" << endl;
  }
  
  void Con::open(string key, int num, char mode) {
    
    struct stat statBuf;
    if(stat(root_.c_str(), &statBuf) != 0) {
      if(-1==mkdir(root_.c_str(), 0777)) {
	throw runtime_error("failed to mkdir root");
      }
    }

    string dirname = (boost::format("%s/%s") % root_ % key).str();
    if(stat(dirname.c_str(), &statBuf) != 0) {
      if(-1==mkdir(dirname.c_str(), 0777)) {
	throw runtime_error("failed to mkdir dir");
      }
    }

    string filepath = (boost::format("%s/%04i.json") % dirname % num).str();
    
    if(mode=='w') {
      file_.open(filepath.c_str(), ios::out);
    } else if(mode=='r') {
      file_.open(filepath.c_str(), ios::in);
    }
    if(!file_.is_open()) 
      throw runtime_error("file open failed");
  }
  void Con::close() {
    file_.close();
  }
  bool Con::has(string key, int num) const {
    string filepath = (boost::format("%s/%s/%04i.json") % root_ % key % num).str();
    ifstream ifs(filepath);
    return ifs.is_open();
  }

  void Con::write_l(string key, int num, bool val) {
    this->open(key, num, 'w');
    if(val){
      file_ << "true" << endl;
    } else {
      file_ << "false" << endl;
    }
    this->close();
  }
  void Con::write_i(string key, int num, int val) {
    this->open(key, num, 'w');
    file_ << val << endl;
    this->close();
  }
  void Con::write_f(string key, int num, double val) {
    this->open(key, num, 'w');
    file_ << val << endl;
    this->close();
  }  
  void Con::write_i1(string key, int num, const VectorXi& val) {
    this->open(key, num, 'w');
    write_vector<int>(file_, val);
    this->close();
  }
  void Con::write_i2(string key, int num, const MatrixXi& val) {
    this->open(key, num, 'w');
    write_matrix<int>(file_, val);
    this->close();
  }
  void Con::write_f1(string key, int num, const VectorXd& val) {
    this->open(key, num, 'w');
    write_vector<double>(file_, val);
    this->close();
  }
  void Con::write_f2(string key, int num, const MatrixXd& val) {
    this->open(key, num, 'w');
    write_matrix<double>(file_, val);
    this->close();
  }
  
}
