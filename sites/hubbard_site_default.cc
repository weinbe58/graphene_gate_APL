#include <dmtk/dmtk.h>

using namespace dmtk;

//////////////////////////////////////////////////////////////////////////////

#if defined(GRANDCANONICAL) 
#warning USING GRANDCANONICAL
template<class T>
void build_hubbard_site(Block<T> &r, double U, double mu, double field = 0, double lambda = 0.0f, T coef = T(0)) 
{
  int dim = 0;
  Basis basis(4);

  QN qn(0);
  basis(0) = State(0,qn);
  qn["N"] = 2;
  basis(1) = State(1,qn); //n = 2
  qn["N"] = 1;
  qn["Sz"] = -1;
  basis(2) = State(2,qn); //n = 1
  qn["Sz"] = 1;
  basis(3) = State(3,qn); // n = 1
  dim = 4;
  basis.reorder();
  for(int i = 0; i < 4; i++) cout << "BASIS " << i << " " << basis(i).qn() << endl;
  r.clear();
  r.resize(basis);
  r.set_lattice(Lattice(1,dmtk::OBC));

  Matrix<T> h(dim,dim);
  Matrix<T> sz(dim,dim), sp(dim,dim), sm(dim,dim), s2(dim,dim);
  Matrix<T> cup(dim,dim), cdup(dim,dim);
  Matrix<T> cdn(dim,dim), cddn(dim,dim);
  Matrix<T> p(dim,dim), pd(dim,dim);
  Matrix<T> n(dim,dim), nup(dim,dim), ndn(dim,dim);
  Matrix<T> dd(dim,dim), dh(dim,dim);

// tables for Hamiltonian

  p(1,0) = 1.;
  pd(0,1) = 1.;

  sz(2,2) = -0.5;
  sz(3,3) = 0.5;

  s2(2,2) = s2(3,3) = 0.75;

  n(2,2) = n(3,3) = 1;
  n(1,1) = 2;

  sm(3,2) = 1;
  sp(2,3) = 1;

  cup(3,0) = 1;
  cup(1,2) = 1;
  cdup(0,3) = 1;
  cdup(2,1) = 1;

  cdn(2,0) = 1;
  cdn(1,3) = -1;
  cddn(0,2) = 1;
  cddn(3,1) = -1;

  dd(1,1) = 1.;
  dh(0,0) = 1.;

  h(2,2) = -mu+0.5*field+lambda*0.75;
  h(3,3) = -mu-0.5*field+lambda*0.75;
  h(1,1) = U - 2.0*mu;

  nup(3,3) = 1;
  nup(1,1) = 1;
  ndn(2,2) = 1;
  ndn(1,1) = 1;

#ifdef USE_CORRELATED_HOPPING
  Matrix<T> cupndn(dim,dim), cdupndn(dim,dim);
  Matrix<T> cdnnup(dim,dim), cddnnup(dim,dim);  
  cupndn = product(cup,ndn);
  cdupndn = product(cdup,ndn);
  cdnnup = product(cdn,nup);
  cddnnup = product(cddn,nup);
#endif

  h += coef*p + dmtk::conj(coef)*pd;

  r.push_back(Sz<T>(sz,r.basis(),0));
  r.push_back(Splus<T>(sp,r.basis(),0));
  r.push_back(Sminus<T>(sm,r.basis(),0));
  r.push_back(Cdup<T>(cdup,r.basis(),0));
  r.push_back(Cddn<T>(cddn,r.basis(),0));
  r.push_back(Cup<T>(cup,r.basis(),0));
  r.push_back(Cdn<T>(cdn,r.basis(),0));
  r.push_back(N<T>(n,r.basis(),0));
  r.push_back(H<T>(h,r.basis(),0));
  r.push_back(Pair<T>(p,r.basis(),0));
  r.push_back(Paird<T>(pd,r.basis(),0));
  r.push_back(Double<T>(dd,r.basis(),0));
  r.push_back(Double<T>(dh,r.basis(),0).set_name("Hole"));
  r.push_back(Nup<T>(nup,r.basis(),0));
  r.push_back(Ndn<T>(ndn,r.basis(),0));
  r.push_back(Sz<T>(s2,r.basis(),0).set_name("S2"));
#ifdef USE_CORRELATED_HOPPING
  r.push_back(Cdup<T>(cdupndn,r.basis(),0).set_name("CdupNdn"));
  r.push_back(Cddn<T>(cddnnup,r.basis(),0).set_name("CddnNup"));
  r.push_back(Cup<T>(cupndn,r.basis(),0).set_name("CupNdn"));
  r.push_back(Cdn<T>(cdnnup,r.basis(),0).set_name("CdnNup"));
#endif

  typename BMatrix<T>::iterator iter;
  Paird<T> op(pd,r.basis(),0);
  for(iter = op.begin(); iter != op.end(); iter++){
    cout << (*iter).qn().n() << " " << (*iter).qn().sz() << endl;
    cout << (*iter).row_range().begin() << " " << (*iter).row_range().end() << " " << (*iter).row_range().size() << endl;
    cout << (*iter).col_range().begin() << " " << (*iter).col_range().end() << " " << (*iter).col_range().size() << endl;
   cout << "--------------------------\n";
   for(int i = 0; i < (*iter).row_range().size(); i++)
     for(int j = 0; j < (*iter).col_range().size(); j++)
       cout << i << " " << j << " " << (*iter)(i,j) << endl;
  }

  cout << "==================================================\n";
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) cout << i << " " << j << " " << op(i,j) << endl;

}
#elif defined(GRANDCANONICAL_SZ) 
#warning USING GRANDCANONICAL_SZ
template<class T>
void build_hubbard_site(Block<T> &r, double U, double mu, double field = 0, double lambda = 0.0f, T coef = T(0)) 
{
  int dim = 0;
  Basis basis(4);

  QN qn(0);
  basis(0) = State(0,qn);
  qn["N"] = 1;
  qn["Sz"] = -1;
  basis(2) = State(1,qn); //n = 1
  qn["Sz"] = 1;
  basis(3) = State(2,qn); // n = 1
  qn["N"] = 2;
  qn["Sz"] = 0;
  basis(1) = State(3,qn); //n = 2
  dim = 4;
  basis.reorder();
  for(int i = 0; i < 4; i++) cout << "BASIS " << i << " " << basis(i).qn() << endl;
  r.clear();
  r.resize(basis);
  r.set_lattice(Lattice(1,dmtk::OBC));

  Matrix<T> h(dim,dim);
  Matrix<T> sz(dim,dim), sp(dim,dim), sm(dim,dim), s2(dim,dim);
  Matrix<T> cup(dim,dim), cdup(dim,dim);
  Matrix<T> cdn(dim,dim), cddn(dim,dim);
  Matrix<T> p(dim,dim), pd(dim,dim);
  Matrix<T> n(dim,dim), nup(dim,dim), ndn(dim,dim);
  Matrix<T> dd(dim,dim), dh(dim,dim);

// tables for Hamiltonian

  p(3,0) = 1.;
  pd(0,3) = 1.;

  sz(1,1) = -0.5;
  sz(2,2) = 0.5;

  s2(1,1) = s2(2,2) = 0.75;

  n(1,1) = n(2,2) = 1;
  n(3,3) = 2;

  sm(2,1) = 1;
  sp(1,2) = 1;

  cup(2,0) = 1;
  cup(3,1) = 1;
  cdup(0,2) = 1;
  cdup(1,3) = 1;

  cdn(1,0) = 1;
  cdn(3,2) = -1;
  cddn(0,1) = 1;
  cddn(2,3) = -1;

  dd(3,3) = 1.;
  dh(0,0) = 1.;

  h(1,1) = -mu+0.5*field+lambda*0.75;
  h(2,2) = -mu-0.5*field+lambda*0.75;
  h(3,3) = U - 2.0*mu;

  nup(3,3) = 1;
  nup(2,2) = 1;
  ndn(3,3) = 1;
  ndn(1,1) = 1;

#ifdef USE_CORRELATED_HOPPING
  Matrix<T> cupndn(dim,dim), cdupndn(dim,dim);
  Matrix<T> cdnnup(dim,dim), cddnnup(dim,dim);
  cupndn = product(cup,ndn);
  cdupndn = product(cdup,ndn);
  cdnnup = product(cdn,nup);
  cddnnup = product(cddn,nup);
#endif

  h += coef*p + dmtk::conj(coef)*pd;

  r.push_back(Sz<T>(sz,r.basis(),0));
  r.push_back(Splus<T>(sp,r.basis(),0));
  r.push_back(Sminus<T>(sm,r.basis(),0));
  r.push_back(Cdup<T>(cdup,r.basis(),0));
  r.push_back(Cddn<T>(cddn,r.basis(),0));
  r.push_back(Cup<T>(cup,r.basis(),0));
  r.push_back(Cdn<T>(cdn,r.basis(),0));
  r.push_back(N<T>(n,r.basis(),0));
  r.push_back(H<T>(h,r.basis(),0));
  r.push_back(Pair<T>(p,r.basis(),0));
  r.push_back(Paird<T>(pd,r.basis(),0));
  r.push_back(Double<T>(dd,r.basis(),0));
  r.push_back(Double<T>(dh,r.basis(),0).set_name("Hole"));
  r.push_back(Nup<T>(nup,r.basis(),0));
  r.push_back(Ndn<T>(ndn,r.basis(),0));
  r.push_back(Sz<T>(s2,r.basis(),0).set_name("S2"));

#ifdef USE_CORRELATED_HOPPING
  r.push_back(Cdup<T>(cdupndn,r.basis(),0).set_name("CdupNdn"));
  r.push_back(Cddn<T>(cddnnup,r.basis(),0).set_name("CddnNup"));
  r.push_back(Cup<T>(cupndn,r.basis(),0).set_name("CupNdn"));
  r.push_back(Cdn<T>(cdnnup,r.basis(),0).set_name("CdnNup"));
#endif

  typename BMatrix<T>::iterator iter;
  Paird<T> op(pd,r.basis(),0);
  for(iter = op.begin(); iter != op.end(); iter++){
    cout << (*iter).qn().n() << " " << (*iter).qn().sz() << endl;
    cout << (*iter).row_range().begin() << " " << (*iter).row_range().end() << " " << (*iter).row_range().size() << endl;
    cout << (*iter).col_range().begin() << " " << (*iter).col_range().end() << " " << (*iter).col_range().size() << endl;
   cout << "--------------------------\n";
   for(int i = 0; i < (*iter).row_range().size(); i++)
     for(int j = 0; j < (*iter).col_range().size(); j++)
       cout << i << " " << j << " " << (*iter)(i,j) << endl;
  }

  cout << "==================================================\n";
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++) cout << i << " " << j << " " << op(i,j) << endl;

}
#elif defined(GRANDCANONICAL_N) 
#warning USING GRANDCANONICAL_N
template<class T>
void build_hubbard_site(Block<T> &r, double U, double mu, double field = 0, double lambda = 0.0f, T coef = T(0)) 
{
  int dim = 0;
  Basis basis(4);

  basis(0) = State(0,QN(1,-1));
  basis(1) = State(1,QN(0,0));
  basis(2) = State(2,QN(2,0)); 
  basis(3) = State(3,QN(1,1)); 
  dim = 4;
  basis.reorder();
  r.clear();
  r.resize(basis);
  r.set_lattice(Lattice(1,dmtk::OBC));

  Matrix<T> h(dim,dim);
  Matrix<T> sz(dim,dim), sp(dim,dim), sm(dim,dim), s2(dim,dim);
  Matrix<T> cup(dim,dim), cdup(dim,dim);
  Matrix<T> cdn(dim,dim), cddn(dim,dim);
  Matrix<T> p(dim,dim), pd(dim,dim);
  Matrix<T> n(dim,dim), nup(dim,dim), ndn(dim,dim);
  Matrix<T> dd(dim,dim), nh(dim,dim);

// tables for Hamiltonian

  p(2,1) = 1.;
  pd(1,2) = 1.;

  sz(0,0) = -0.5;
  sz(3,3) = 0.5;

  s2(0,0) = s2(3,3) = 0.75;

  n(0,0) = n(3,3) = 1;
  n(2,2) = 2;

  sm(3,0) = 1;
  sp(0,3) = 1;

  cup(2,0) = 1;
  cup(3,1) = 1;
  cdup(0,2) = 1;
  cdup(1,3) = 1;

  cdn(2,3) = -1;
  cdn(0,1) = 1;
  cddn(3,2) = -1;
  cddn(1,0) = 1;

  dd(2,2) = 1.;

  h(0,0) = -mu+0.5*field+lambda*0.75;
  h(3,3) = -mu-0.5*field+lambda*0.75;
  h(2,2) = U - 2.0*mu;

  nup(3,3) = 1;
  nup(2,2) = 1;
  ndn(2,2) = 1;
  ndn(0,0) = 1;

  nh = T(0);
  nh(0,0) = 1;

  h += coef*p + dmtk::conj(coef)*pd;
#ifdef USE_CORRELATED_HOPPING
  Matrix<T> cupndn(dim,dim), cdupndn(dim,dim);
  Matrix<T> cdnnup(dim,dim), cddnnup(dim,dim);
  cupndn = product(cup,ndn);
  cdupndn = product(cdup,ndn);
  cdnnup = product(cdn,nup);
  cddnnup = product(cddn,nup);
#endif

  r.push_back(Sz<T>(sz,r.basis(),0));
  r.push_back(Splus<T>(sp,r.basis(),0));
  r.push_back(Sminus<T>(sm,r.basis(),0));
  r.push_back(Cdup<T>(cdup,r.basis(),0));
  r.push_back(Cddn<T>(cddn,r.basis(),0));
  r.push_back(Cup<T>(cup,r.basis(),0));
  r.push_back(Cdn<T>(cdn,r.basis(),0));
  r.push_back(N<T>(n,r.basis(),0));
  r.push_back(N<T>(nh,r.basis(),0).set_name("Nh"));
  r.push_back(N<T>(n,r.basis(),0));
  r.push_back(H<T>(h,r.basis(),0));
  r.push_back(Pair<T>(p,r.basis(),0));
  r.push_back(Paird<T>(pd,r.basis(),0));
  r.push_back(Double<T>(dd,r.basis(),0));
  r.push_back(Nup<T>(nup,r.basis(),0));
  r.push_back(Ndn<T>(ndn,r.basis(),0));
  r.push_back(Sz<T>(s2,r.basis(),0).set_name("S2"));

#ifdef USE_CORRELATED_HOPPING
  Matrix<T> cupndn(dim,dim), cdupndn(dim,dim);
  Matrix<T> cdnnup(dim,dim), cddnnup(dim,dim);
  cupndn = product(cup,ndn);
  cdupndn = product(cdup,ndn);
  cdnnup = product(cdn,nup);
  cddnnup = product(cddn,nup);
  r.push_back(Cdup<T>(cdupndn,r.basis(),0).set_name("CdupNdn"));
  r.push_back(Cddn<T>(cddnnup,r.basis(),0).set_name("CddnNup"));
  r.push_back(Cup<T>(cupndn,r.basis(),0).set_name("CupNdn"));
  r.push_back(Cdn<T>(cdnnup,r.basis(),0).set_name("CdnNup"));
#endif
/*
  typename BMatrix<T>::iterator iter;
  Cup<T> op(cup,r.basis(),0);
  for(iter = op.begin(); iter != op.end(); iter++){
    cout << (*iter).qn().n() << " " << (*iter).qn().sz() << endl;
    cout << (*iter).row_range().begin() << " " << (*iter).row_range().end() << " " << (*iter).row_range().size() << endl;
    cout << (*iter).col_range().begin() << " " << (*iter).col_range().end() << " " << (*iter).col_range().size() << endl;
   for(int i = 0; i < (*iter).row_range().size(); i++)
     for(int j = 0; j < (*iter).col_range().size(); j++)
       cout << i << " " << j << " " << (*iter)(i,j) << endl;
  }
*/

}
#else // USE_PAIR

template<class T>
void build_hubbard_site(Block<T> &r, double U, double mu, double field = 0, double lambda = 0, T coef = T(0))
{   
  int dim = 0;
  Basis basis(4);
 
  QN qn(0); 
  basis(0) = State(0,qn);
  qn["N"]=1;
  qn["Sz"]=-1;
  basis(1) = State(1,qn);
  qn["N"]=1;
  qn["Sz"]=1;
  basis(2) = State(2,qn);
  qn["N"]=2;
  qn["Sz"]=0;
  basis(3) = State(3,qn);
  dim = 4;
  basis.reorder();
  r.clear();
  r.resize(basis);
  r.set_lattice(Lattice(1,dmtk::OBC));

  Matrix<T> h(dim,dim);
  Matrix<T> sz(dim,dim), sp(dim,dim), sm(dim,dim), s2(dim,dim);
  Matrix<T> cup(dim,dim), cdup(dim,dim);
  Matrix<T> cdn(dim,dim), cddn(dim,dim);
  Matrix<T> p(dim,dim), pd(dim,dim);
  Matrix<T> n(dim,dim), nup(dim,dim), ndn(dim,dim);
  Matrix<T> dd(dim,dim), dh(dim,dim);
  Matrix<T> ns(dim,dim), ne(dim,dim);

// tables for Hamiltonian

  p(3,0) = 1.;
  pd(0,3) = 1.;

  sz(1,1) = -0.5;
  sz(2,2) = 0.5;
 
  s2(1,1) = s2(2,2) = 0.75;

  n(1,1) = n(2,2) = 1;
  n(3,3) = 2;

  ndn(1,1) = 1;
  ndn(3,3) = 1;
  nup(2,2) = 1;
  nup(3,3) = 1;

  ns(0,0) = 0;
  ns(1,1) = 1;
  ns(2,2) = 1;
  ns(3,3) = 0;

  ne(0,0) = 1;

  sm(2,1) = 1;
  sp(1,2) = 1;

  cup(2,0) = 1;
  cup(3,1) = 1;
  cdup(0,2) = 1;
  cdup(1,3) = 1;

  cdn(1,0) = 1;
  cdn(3,2) = -1;
  cddn(0,1) = 1;
  cddn(2,3) = -1;

  dd(3,3) = 1.;
  dh(0,0) = 1.;

  h(1,1) = -mu-0.5*field+lambda*0.75;
  h(2,2) = -mu+0.5*field+lambda*0.75;
  h(3,3) = U - 2.0*mu;
  h += coef*p + dmtk::conj(coef)*pd;

#ifdef USE_CORRELATED_HOPPING
  Matrix<T> cupndn(dim,dim), cdupndn(dim,dim);
  Matrix<T> cdnnup(dim,dim), cddnnup(dim,dim);
  cupndn = product(cup,ndn);
  cdupndn = product(cdup,ndn);
  cdnnup = product(cdn,nup);
  cddnnup = product(cddn,nup);
#endif

  r.push_back(Sz<T>(sz,r.basis(),0));
  r.push_back(Splus<T>(sp,r.basis(),0));
  r.push_back(Sminus<T>(sm,r.basis(),0));
  r.push_back(Cdup<T>(cdup,r.basis(),0));
  r.push_back(Cddn<T>(cddn,r.basis(),0));
  r.push_back(Cup<T>(cup,r.basis(),0));
  r.push_back(Cdn<T>(cdn,r.basis(),0));
  r.push_back(Double<T>(dd,r.basis(),0));
  r.push_back(Double<T>(dh,r.basis(),0).set_name("Hole"));
  r.push_back(Nup<T>(nup,r.basis(),0));
  r.push_back(Ndn<T>(ndn,r.basis(),0));
  r.push_back(N<T>(n,r.basis(),0));
  r.push_back(N<T>(ns,r.basis(),0).set_name("Nsingle"));
  r.push_back(N<T>(ne,r.basis(),0).set_name("Nempty"));
  r.push_back(H<T>(h,r.basis(),0));
  r.push_back(Pair<T>(p,r.basis(),0));
  r.push_back(Paird<T>(pd,r.basis(),0));
  r.push_back(Sz<T>(s2,r.basis(),0).set_name("S2"));
#ifdef USE_CORRELATED_HOPPING
  r.push_back(Cdup<T>(cdupndn,r.basis(),0).set_name("CdupNdn"));
  r.push_back(Cddn<T>(cddnnup,r.basis(),0).set_name("CddnNup"));
  r.push_back(Cup<T>(cupndn,r.basis(),0).set_name("CupNdn"));
  r.push_back(Cdn<T>(cdnnup,r.basis(),0).set_name("CdnNup"));
#endif
}

#endif // GRANDCANONICAL


