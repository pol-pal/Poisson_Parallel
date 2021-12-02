template<typename Number>
class Problem
{
    Number eta;
public:
  typedef Number value_type;

  //! Constructor without arg sets nonlinear term to zero
  Problem () : eta(0.0) {}

 //! Constructor takes eta parameter
  Problem (const Number& eta_) : eta(eta_) {}

 
//! right hand side
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    Number s=1.0; 
	for (std::size_t i=0; i<global.size(); i++){
          s*=sin(M_PI*global[i]); 
	 }
         return M_PI*M_PI*s*x.size();
  }

  //! boundary condition type function (true = Dirichlet)
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return true;
  }

  //! Dirichlet extension
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    Number s=1.0;
    for (std::size_t i=0; i<global.size(); i++){ 
       s*=sin(M_PI*global[i]);
    }
    return s;
  }


};


