#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

template<typename F, typename FiniteElementMap>
class PoissonPk :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,

  public Dune::PDELab::
    NumericalJacobianVolume<PoissonPk<F,FiniteElementMap> >,
  public Dune::PDELab::
    NumericalJacobianApplyVolume<PoissonPk<F,FiniteElementMap> >
{
private:
  // define useful types
  typedef typename FiniteElementMap::Traits::FiniteElementType
     FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType
     LocalBasisType;
  typedef typename LocalBasisType::Traits::DomainType
     DomainType;
  typedef typename LocalBasisType::Traits::RangeFieldType
     RF;
  typedef typename LocalBasisType::Traits::RangeType
     RangeType;
  typedef typename LocalBasisType::Traits::JacobianType
     JacobianType;
   
  Dune::PDELab::LocalBasisCache<LocalBasisType> cache;
  
// data members
  enum {dim=LocalBasisType::Traits::dimDomain};
  const F f; 

public:
  // define flags controlling global assembler
  enum { doPatternVolume = true }; //determine the sparsity of the matrix A
  enum { doAlphaVolume = true }; //volume integral  for solution u_h
 // enum { doLambdaVolume = true };//volume integral for RHS

  // Constructor precomputes element independent data
  PoissonPk (const F& f_, const FiniteElementType& fel)
    : f(f_)
  {
    
  }


// jacobian of volume term
  //compute the element contributions to the stiffness matrix
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu,
                        const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    
    size_t n=lfsu.size(); //number of basis functions

    // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > grad(lfsu.size());
  
   // select quadrature rule
    auto geo = eg.geometry();
    const int order =2*lfsv.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);
    
   

 // loop over quadrature points
    for (const auto& ip : rule){

   //compute common factor for every quadrature point factor=w_d*|detB_t|
      RF factor = ip.weight()*geo.integrationElement(ip.position()); 

    // evaluate basis functions
        auto& phihat = cache.evaluateFunction(ip.position(),
                             lfsu.finiteElement().localBasis());

      // evaluate gradient of shape functions
    auto& gradhat = cache.evaluateJacobian(ip.position(),
                             lfsu.finiteElement().localBasis());

    // transform gradients of shape functions to real element
        const auto S = geo.jacobianInverseTransposed(ip.position()); 
        for (size_t i=0; i<lfsu.size(); i++)
           S.mv(gradhat[i][0],grad[i]);
   
    // store in result
    for (int i=0; i<n; i++)
         for (int j=0; j<n; j++)
            mat.accumulate(lfsu,i,lfsu,j,(grad[i]*grad[j]*factor));//accumulate result for local stiffness matrix  A_T=grad^T*grad*factor
  }
    
}

// volume integral depending on test and ansatz functions
  //element local computations for matrix free evaluations of a(u_h,phi_i) for all test functions phi_i
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu,
                     const X& x, const LFSV& lfsv,
                     R& r) const
{ 

   typename F::Traits::RangeType fval;

   size_t n=lfsu.size();

   // Initialize vectors outside for loop
        std::vector<Dune::FieldVector<RF,dim> > grad(lfsu.size());
        Dune::FieldVector<RF,dim> graduh(0.0);

// select quadrature rule
    auto geo = eg.geometry();
    const int order = 2*lfsv.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

 // loop over quadrature points
    for (const auto& ip : rule){  

      // evaluate basis functions
        auto& phihat = cache.evaluateFunction(ip.position(),
                             lfsu.finiteElement().localBasis());

      // evaluate gradient of shape functions
      auto& gradhat = cache.evaluateJacobian(ip.position(),
                             lfsu.finiteElement().localBasis());
 
     //evaluate RHS at quadratic point
         f.evaluate(eg.entity(),ip.position(),fval);

     //compute common factor1= f*w_d*|detB_T|
         RF factor1=fval*ip.weight()*eg.geometry().integrationElement(ip.position());

     // transform gradients of shape functions to real element
        const auto S = geo.jacobianInverseTransposed(ip.position());

    //    auto grad = Dune::PDELab::makeJacobianContainer(lfsu);
        for (size_t i=0; i<lfsu.size(); i++)
          S.mv(gradhat[i][0],grad[i]);

     // compute gradient of u
        graduh = 0.0;
        for (size_t i=0; i<lfsu.size(); i++)
          graduh.axpy(x(lfsu,i),grad[i]);

    // integrate (grad u)*grad phi_i
       auto factor2 = ip.weight()*
          geo.integrationElement(ip.position());

       for (size_t i=0; i<lfsu.size(); i++)
          r.accumulate(lfsu,i,(graduh*grad[i])*factor2-(factor1*phihat[i]));
    }
}

  //! apply local jacobian of the volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }

  //! apply local jacobian of the volume term
  template<typename EG, typename LFSU, typename X,
           typename LFSV, typename R>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu,
                              const X& x, const X& z, const LFSV& lfsv,
                              R& r) const
  {
    alpha_volume(eg,lfsu,z,lfsv,r);
  }
};
