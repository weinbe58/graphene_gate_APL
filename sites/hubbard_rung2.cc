#include <dmtk/dmtk.h>

namespace dmtk{

    
template<class T>
class CddnCup: public dmtk::BasicOp<T>
    {
    private:
        typedef BasicOp<T> _O;
        virtual void init()
        { _O::dqn.sz() = -2; _O::dqn.n() = 0; _O::_is_diagonal = false; }
    public:
        CddnCup(): BasicOp<T>("CddnCup",0)
        { init(); }
        CddnCup(int site, int internal_site = 0): BasicOp<T>("CddnCup",site,internal_site)  { init(); }
        CddnCup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("CddnCup",b,site,internal_site)
            { init(); _O::init_blocks(); }
        CddnCup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("CddnCup",m,b,site,internal_site)
            { init(); bmatrix_from_matrix(m); }
            
    };
    
    
    
template<class T>
class CdupCdn: public dmtk::BasicOp<T>
    {
    private:
        typedef BasicOp<T> _O;
        virtual void init()
        { _O::dqn.sz() = +2; _O::dqn.n() = 0; _O::_is_diagonal = false; }
    public:
        CdupCdn(): BasicOp<T>("CdupCdn",0)
        { init(); }
        CdupCdn(int site, int internal_site = 0): BasicOp<T>("CdupCdn",site,internal_site)  { init(); }
        CdupCdn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("CdupCdn",b,site,internal_site)
        { init(); _O::init_blocks(); }
        CdupCdn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("CdupCdn",m,b,site,internal_site)
        { init(); bmatrix_from_matrix(m); }
        
    };
    
    
template<class T>
class CupCddn: public dmtk::BasicOp<T>
    {
    private:
        typedef BasicOp<T> _O;
        virtual void init()
        { _O::dqn.sz() = -2; _O::dqn.n() = 0; _O::_is_diagonal = false; }
    public:
        CupCddn(): BasicOp<T>("CupCddn",0)
        { init(); }
        CupCddn(int site, int internal_site = 0): BasicOp<T>("CupCddn",site,internal_site)  { init(); }
        CupCddn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("CupCddn",b,site,internal_site)
        { init(); _O::init_blocks(); }
        CupCddn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("CupCddn",m,b,site,internal_site)
        { init(); bmatrix_from_matrix(m); }
        
    };
    
template<class T>
class CdnCdup: public dmtk::BasicOp<T>
    {
    private:
        typedef BasicOp<T> _O;
        virtual void init()
        { _O::dqn.sz() = +2; _O::dqn.n() = 0; _O::_is_diagonal = false; }
    public:
        CdnCdup(): BasicOp<T>("CdnCdup",0)
        { init(); }
        CdnCdup(int site, int internal_site = 0): BasicOp<T>("CdnCdup",site,internal_site)  { init(); }
        CdnCdup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("CdnCdup",b,site,internal_site)
        { init(); _O::init_blocks(); }
        CdnCdup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("CdnCdup",m,b,site,internal_site)
        { init(); bmatrix_from_matrix(m); }
        
    };

template<class T>
class CdupCup: public dmtk::BasicOp<T>
    {
    private:
        typedef BasicOp<T> _O;
        virtual void init()
        { _O::dqn.sz() = +0; _O::dqn.n() = 0; _O::_is_diagonal = false; }
    public:
        CdupCup(): BasicOp<T>("CdupCup",0)
        { init(); }
        CdupCup(int site, int internal_site = 0): BasicOp<T>("CdupCup",site,internal_site)  { init(); }
        CdupCup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("CdupCup",b,site,internal_site)
        { init(); _O::init_blocks(); }
        CdupCup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("CdupCup",m,b,site,internal_site)
        { init(); bmatrix_from_matrix(m); }

    };

template<class T>
class CddnCdn: public dmtk::BasicOp<T>
    {
    private:
        typedef BasicOp<T> _O;
        virtual void init()
        { _O::dqn.sz() = +0; _O::dqn.n() = 0; _O::_is_diagonal = false; }
    public:
        CddnCdn(): BasicOp<T>("CddnCdn",0)
        { init(); }
        CddnCdn(int site, int internal_site = 0): BasicOp<T>("CddnCdn",site,internal_site)  { init(); }
        CddnCdn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("CddnCdn",b,site,internal_site)
        { init(); _O::init_blocks(); }
        CddnCdn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("CddnCdn",m,b,site,internal_site)
        { init(); bmatrix_from_matrix(m); }

    };

    

} // namespace dmtk

using namespace dmtk;

template<class T>
void build_2sites(Block<T> &r, double U, double t_rung = 1.0, double xmu = 0.0f, T x = T(0)) 
{
  Block<T> aux;
  build_hubbard_site(aux, U, xmu);
  Basis aux_basis(aux.basis(),aux.basis());
  aux_basis.reorder();

  r.clear();
  r.resize(aux_basis);
  r.set_lattice(Lattice(1,OBC));
  r.set_orbitals(2);

  BMatrix<T> aux_rho(aux_basis);
  PackedBasis::const_iterator pbiter;

  for(pbiter = aux_basis.subspaces().begin(); pbiter != aux_basis.subspaces().end(); pbiter++){
    SubMatrix<T> block(pbiter->qn(),*pbiter,*pbiter);
    block = I<T>();
    block.qn() = (*pbiter).qn();
    aux_rho.push_back(block);
  }

  N<T> n0(aux_basis,0);
  new_operator(n0, *aux("N",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(n0);
  Sz<T> sz0(aux_basis,0);
  new_operator(sz0, *aux("Sz",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(sz0);
  Splus<T> sp0(aux_basis,0);
  new_operator(sp0, *aux("S+",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(sp0);
  Sminus<T> sm0(aux_basis,0);
  new_operator(sm0, *aux("S-",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(sm0);

  Cup<T> cup0(aux_basis,0);
  new_operator(cup0, *aux("Cup",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(cup0);
  Cdup<T> cdup0(aux_basis,0);
  new_operator(cdup0, *aux("Cdup",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(cdup0);
  Cdn<T> cdn0(aux_basis,0);
  new_operator(cdn0, *aux("Cdn",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(cdn0);
  Cddn<T> cddn0(aux_basis,0);
  new_operator(cddn0, *aux("Cddn",0), aux_rho, aux_basis, LEFT, T(1), false);
  r.push_back(cddn0);

  N<T> n1(aux_basis,0,1);
  new_operator(n1, *aux("N",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(n1);
  Sz<T> sz1(aux_basis,0,1);
  new_operator(sz1, *aux("Sz",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(sz1);
  Splus<T> sp1(aux_basis,0,1);
  new_operator(sp1, *aux("S+",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(sp1);
  Sminus<T> sm1(aux_basis,0,1);
  new_operator(sm1, *aux("S-",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(sm1);

  Cup<T> cup1(aux_basis,0,1);
  new_operator(cup1, *aux("Cup",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(cup1);
  Cdup<T> cdup1(aux_basis,0,1);
  new_operator(cdup1, *aux("Cdup",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(cdup1);
  Cdn<T> cdn1(aux_basis,0,1);
  new_operator(cdn1, *aux("Cdn",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(cdn1);
  Cddn<T> cddn1(aux_basis,0,1);
  new_operator(cddn1, *aux("Cddn",0), aux_rho, aux_basis, RIGHT, T(1), false);
  r.push_back(cddn1);


  // Added
    
    CddnCup<T>  cddncup(aux_basis, 0);
    new_operator(cddncup, *aux("Cddn",0), *aux("Cup",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
    r.push_back(cddncup);

    CdupCdn<T>  cdupcdn(aux_basis, 0);
    new_operator(cdupcdn, *aux("Cdup",0), *aux("Cdn",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
    r.push_back(cdupcdn);

    CupCddn<T>  cupcddn(aux_basis, 0);
    new_operator(cupcddn, *aux("Cup",0), *aux("Cddn",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
    r.push_back(cupcddn);

    CdnCdup<T>  cdncdup(aux_basis, 0);
    new_operator(cdncdup, *aux("Cdn",0), *aux("Cdup",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
    r.push_back(cdncdup);

    CddnCdn<T>  cddncdn(aux_basis, 0);
    new_operator(cddncdn, *aux("Cddn",0), *aux("Cdn",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
    r.push_back(cddncdn);

    CdupCup<T>  cdupcup(aux_basis, 0);
    new_operator(cdupcup, *aux("Cdup",0), *aux("Cup",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
    r.push_back(cdupcup);

 
  Sz<T> szsz(aux_basis, 0);
  szsz.set_name("Sz(0,0)Sz(0,1)");
  new_operator(szsz, *aux("Sz",0), *aux("Sz",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1), false, false);
  r.push_back(szsz);



  typename Block<T>::iterator iter;
  for(iter = aux.begin(); iter != aux.end(); iter++)
    cout << "OPERATOR " << (*iter).name() << " " << (*iter).site() << " " << (*iter).internal_site() << endl;



/*

  H<T> h(aux_basis,0);
#ifdef USE_SITE_PAIR
  new_operator(h, *aux("Paird",0), aux_rho, aux_basis, LEFT, T(x), false);
  new_operator(h, *aux("Pair",0), aux_rho, aux_basis, LEFT, dmtk::conj(x), false);
#endif
  new_operator(h, *aux("Paird",0), aux_rho, aux_basis, RIGHT, T(x), false);
  new_operator(h, *aux("Pair",0), aux_rho, aux_basis, RIGHT, dmtk::conj(x), false);
 
  new_operator(h, *aux(H<T>()), aux_rho, aux_basis, LEFT, T(1), false);
  new_operator(h, *aux(H<T>()), aux_rho, aux_basis, RIGHT, T(1), false);
  new_operator(h, *aux("Cdup",0), *aux("Cup",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(-t_rung), true, false);
  new_operator(h, *aux("Cddn",0), *aux("Cdn",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(-t_rung), true, false);

#ifndef USE_SITE_PAIR
  new_operator(h, *aux("Cup",0), *aux("Cdn",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(-1./DMTK_SQRT2)*x, true, false);
  new_operator(h, *aux("Cdn",0), *aux("Cup",0), BLOCK1, BLOCK2, aux_rho, aux_basis, T(1./DMTK_SQRT2)*x, true, false);
#endif

  r.push_back(h);

}

*/
}
