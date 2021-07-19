#define USE_VDISK
// #define WITH_NR
#define WIDTH 25
// #define USE_ANDERSON
#include <fstream>
#include <string>
#include <stdexcept>
#include <dmtk/dmtk.h>
#include <dmtk/site.h>
#include "sites/hubbard_site_default.cc"
#include "sites/spin_site.cc"



using namespace dmtk;
// int impurity_site = 0;
// int lc = 0;
Block<double> impurity_block_left,impurity_block_right;


void
read_mapping(string, int,Vector<double>&, Vector<double>&);

template<class T>
Hami<T>
impurity(int lx, Vector<double>& b, Vector<double>& a, double U);

int main()
{
  double U, V=0, Vg=0, mu_chain, t, mu_gr, gap;
  int lx, lc, radgr, ix;
  int nexp;
  int ntrunc, m;
  int niter;
  size_t dir = RIGHT2LEFT;
  int start = 1;
  int custom;
  int nsteps;
  int nstates[2][20];
  int iop;
  string mapping_file;

                                                                             
  cout << "Input chain length: ";
  cin >> lc;
  cout << lc << endl;

  cout << "Input graphene radius: ";
  cin >> radgr;
  cout << radgr << endl;

  lc = max(0,lc);
  radgr = max(0,radgr);

  lx = 2*radgr+lc+2; // two graphene flakes with the chain bridge and two extra sites for the impurities

  cout << "Input chain chemical potential: ";
  cin >> mu_chain;
  cout << mu_chain << endl;

  cout << "Input chain hopping: ";
  cin >> t;
  cout << t << endl;

  cout << "Input chain gap: ";
  cin >> gap;
  cout << gap << endl;

  cout << "Input Graphene chemical potential: ";
  cin >> mu_gr;
  cout << mu_gr << endl;

#ifdef USE_ANDERSON
  cout << "Input U: ";
  cin >> U;
  cout << U << endl;

  cout << "Input Vg: ";
  cin >> Vg;
  cout << Vg << endl;

  cout << "Input V: ";
  cin >> V;
  cout << V << endl;
#else
  Vg = 0;
  V = 0;
  cout << "Input Jk: ";
  cin >> U;
  cout << U << endl;
#endif

  QN::add_qn_index("N",true);
  QN::add_qn_index("Sz",false);
  QN::set_qn_mask(QN::default_mask());

  int xn,xsz;
  cout << "Input Nf: ";
  cin >> xn;
  cout << xn << endl;

  cout << "mapping file: ";
  cin >> mapping_file;
  cout << mapping_file << endl;

  // cout << "Input Sz: ";
  // cin >> xsz;
  // cout << xsz << endl;



  cout << "Customize dmrg? (1->YES // 0->NO): ";
  cin >> custom;
  cout << custom  << endl;

  if(custom == 0){
    cout << "Number of states: ";
    cin >> ntrunc;
    cout << ntrunc << endl;
                                                                             
    m = ntrunc/6;
    nstates[0][0] = 20; // warmup
    nstates[0][1] = 2*m; nstates[1][1] = 3*m;
    nstates[0][2] = 4*m; nstates[1][2] = 5*m;
    nstates[0][3] = 6*m; nstates[1][3] = 6*m;
    nstates[0][4] = 6*m; nstates[1][4] = 6*m;
    nsteps = 4;
  } else {
    cout << "Number of dmrg iterations:" << endl;
    cin >> nsteps;
                                                                             
    for(int i = 1; i <= nsteps; i++){
      cout << "Iteration " << i << " :" << endl;
      cout << "Number of states L2R:" << endl;
      cin >> nstates[0][i];
      cout << "Number of states R2L:" << endl;
      cin >> nstates[1][i];
    }
  }
  cout << "Start from:" << endl;
  cout << "    0 - Warmup loop" << endl;
  for(int i = 1; i <= nsteps; i++){
    cout << "    " << i << " - Iteration " << i << " : (" << nstates[0][i] << "," << nstates[1][i] << ") states" << endl;
  }
                                                                             
  cout << "    999 - Final Iteration" << endl;
  cin >> niter;
                                                                             
  if(niter >= 1 && niter <= nsteps || niter == 999){
    cout << "Direction (0->LEFT2RIGHT // 1->RIGHT2LEFT):\n";
    cin >> dir;
    cout << "Iteration:\n";
    cin >> start;
  }

  // reading in mapping file and constructing chain SP description

  Vector<double> a_gr(radgr+1,1),b_gr(radgr+1,1),a(lx,1),b(lx,1);

  if(radgr>0){  
    try
    {
      read_mapping(mapping_file,radgr,a_gr,b_gr); //cause an exception to throw
    }
    catch (invalid_argument& e)
    {
      cerr << e.what() << endl;
      return -1;
    }

    b_gr(radgr-1) = t;
  }



  // set up a:
  ix = 0;
  a(ix++) = -Vg;

  for(int i=0;i<radgr;i++){
    a(ix++) = a_gr(i)-mu_gr;
  }

  for(int i=0;i<lc;i++){
    a(ix++) = -mu_chain;
  }

  for(int i=0;i<radgr;i++){
    a(ix++) = a_gr(radgr-i-1)-mu_gr;
  }

  a(ix++) = -Vg;


  // set up b:
  ix = 0;

  b(ix++) = V;

  for(int i=0;i<radgr-1;i++){
    b(ix++) = b_gr(i);
  }

  b(ix++) = t;

  for(int i=0;i<lc-1;i++){
    b(ix++) = t+(gap/4)*(i%2 ? -1 : 1);
  } 

  if(lc > 0){
    b(ix++) = t;    
  }

  for(int i=0;i<radgr-1;i++){
    b(ix++) = b_gr(radgr-i-2);
  }

  b(ix++) = V;






  Hami<double>hami = impurity<double>(lx, a, b, U);
  System<double> S(hami, hami.lattice(), "impurity");

  S.set_lanczos_tolerance(1.e-9);
  S.set_grow_symmetric(false);
  S.set_calc_gap(0);
  S.set_use_hc(true);

  S.set_store_products(false);
  S.set_qn_mask(QN::default_mask());
  S.qnt["N"] = xn+2;
  S.qnt["Sz"] = (xn%2);
                                              
  cout << "**********************************\n";


  niter = 0;                                                                             
  if(niter == 0){
    S.warmup_loop(40);
    niter = 1;
    cout << "**********************************\n";
  }
  if(niter != 999){
    for(int i = niter; i < nsteps; i++){
      cout << "NEW SWEEP " << i << " " << nstates[0][i] << " " << nstates[1][i] << endl;
      S.sweep(nstates[0][i],nstates[1][i],(size_t)dir,start);
      dir = RIGHT2LEFT;
      start = 1;
    }
  }

  S.set_error(1.e-10);
  S.dm_boundary_ops.clear();
  S.sweep(ntrunc,ntrunc,(size_t)dir,start);
  S.sweep(ntrunc,ntrunc,(size_t)dir,start);
                                                                             
  cout << "**********************************\n";
  cout << "Creating measurement operators\n";
  cout << "**********************************\n";
// Operators to measure
  S.corr.set_use_hc(true);
  for(int i=1;i<lx;i++){
    S.corr += Splus<double>(0)*Sminus<double>(i)*0.5;
    S.corr += Sz<double>(0)*Sz<double>(i);
  }

  for(int i=1;i<lx-1;i++){
    S.corr += N<double>(i);
    S.corr += Nup<double>(i);
    S.corr += Ndn<double>(i);
  }

  S.final_sweep(nstates[1][nsteps],(size_t)dir,start,true);

  S.measure_n(1);
  S.measure();

  cout << "**********************************\n";
  cout << " Printing Results of Measurements \n";
  cout << "**********************************\n";



  cout << setprecision(14) << scientific;
  cout << "RESULTS {";
  cout << "\"Nf\":" << xn << ",";
  cout << "\"E\":" << S.energy[0] << ",";
  for(auto corr_iter = S.corr.begin(); corr_iter != S.corr.end(); corr_iter++){
     cout << "\"" << corr_iter->description() << "\":" << corr_iter->value() << ",";
  }
  cout << "}\n";

 
}

void
read_mapping(string mapping_file,int rad,Vector<double>& A, Vector<double>& B){
  ifstream infile(mapping_file.c_str());
  std::string line;
  int index = 0;
  bool finished = false;


  if(!infile.is_open()){
  	cout << "RUNTIME ERROR: failed to open mapping file" << endl;
    exit(-1);
  }

  while (getline(infile, line))
  {
    istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) {
      break;
    } // error

    if(index < rad-1){
	    A(index) = a;
	    B(index) = b;
    }
    else{
	    A(index) = a;
    	finished = true;
    	break;
    }
	index++;
  }

  if(!finished){
  	throw invalid_argument("radgr too large for mapping!");
  }
}



template<class T>
Hami<T>
impurity(int lx, Vector<double>& a, Vector<double>& b, double U)
{
  Hami<T> chain(Lattice(lx,OBC));
  Lattice l = chain.lattice();
  chain.set_use_hc(true);

   
#ifdef USE_ANDERSON
  chain += Double<T>(0)*T(U);
  chain += Double<T>(lx-1)*T(U);
  for (int ix=0; ix < lx-1; ix++) //for b terms
  {
    chain += Cdup<T>(ix)*Cup<T>(ix+1)*T(-b(ix));
    chain += Cddn<T>(ix)*Cdn<T>(ix+1)*T(-b(ix));
  }
  for (int ix=0; ix < lx; ix++) //for a terms
  {
    chain += N<T>(ix)*T(a(ix));
  }
#else
  for (int ix=1; ix < lx-2; ix++) //for b terms
  {
    chain += Cdup<T>(ix)*Cup<T>(ix+1)*T(-b(ix));
    chain += Cddn<T>(ix)*Cdn<T>(ix+1)*T(-b(ix));
  }
  for (int ix=1; ix < lx-1; ix++) //for a terms
  {
    chain += N<T>(ix)*T(a(ix));
  }

  build_spin_site(impurity_block_left,1);
  chain.sites[0] = &impurity_block_left;

  build_spin_site(impurity_block_right,1);
  chain.sites[lx-1] = &impurity_block_right;

  // chain += Sz<T>(0)*1e-7;
  // chain += Sz<T>(lx-1)*1e-7;
  chain += Splus<T>(0)*Sminus<T>(1)*U*0.5;
  chain += Sz<T>(0)*Sz<T>(1)*U;
  chain += Splus<T>(lx-2)*Sminus<T>(lx-1)*U*0.5;
  chain += Sz<T>(lx-2)*Sz<T>(lx-1)*U;

#endif

  build_hubbard_site<T>(chain.site, 0, 0);

  cout << chain.description() << endl;

  return chain;
}


