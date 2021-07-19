#include <dmtk/dmtk.h>

using namespace dmtk;

static dmtk::Vector<Block<double> > site_block;

template<class T>
void build_spin_site(Block<T> &r, int s) 
{
  double S = s*0.5;
  int dim = 2*S+1.001;
  Basis basis(dim);

  QN qn(0);
  
  if(QN::get_qn_index("N") != QN_MAX_SIZE) qn["N"] = 1;

  int n = 0;
  for(int m = -s; m <= s; m+=2){
    qn["Sz"] = m;
    basis(n) = State(n,qn);
    cout << "BASIS " << n << " " << qn << endl;
    n++;
  }
  basis.reorder();

  r.clear();
  r.resize(basis);
  r.set_lattice(Lattice(1,OBC));

  Matrix<T> h(dim,dim);
  Matrix<T> sz(dim,dim), sp(dim,dim), sm(dim,dim);

// tables for Hamiltonian
  if(s == 1){
    Matrix<T> nup(dim,dim);
    Matrix<T> ndn(dim,dim);
    nup = 0.;
    ndn = 0.;
    ndn(0,0) = 1.;
    nup(1,1) = 1.;
    r.push_back(N<T>(nup,r.basis(),0).set_name("Nup"));
    r.push_back(N<T>(ndn,r.basis(),0).set_name("Ndn"));
  }

  n = 0;
  for(int m = -s; m <= s; m+=2){
    double Sz = m*0.5; 
    sz(n,n) = Sz;
    if(n < dim-1){
      sp(n,n+1) = sqrt(S*(S+1)-Sz*(Sz+1));
      sm(n+1,n) = sqrt(S*(S+1)-Sz*(Sz+1));
      cout << "MATRIX " << n << " " << sz(n,n) << " " << sp(n,n+1) << endl;
    }
    n++;
  }

  r.push_back(Sz<T>(sz,r.basis(),0));
  r.push_back(Splus<T>(sp,r.basis(),0));
  r.push_back(Sminus<T>(sm,r.basis(),0));
  r.push_back(H<T>(h,r.basis(),0));
}

