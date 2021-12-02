#include <dune/pdelab/test/l2difference.hh>
#include "problemCube.hh"
//#include "problemSphere.hh"
//#include"problemTetrahedral.hh"
#include <vector>

//#include "problem.hh"
template<typename GV, typename FEM> //class take the mesh, the fem map, and access to the input file
void driver (const GV& gv, const FEM& fem,
             Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF; // type for ccordinates
  typedef double RF;                   // type for computations

//make PDE parameter class
  RF eta = ptree.get("problem.eta",(RF)1.0);
  Problem<RF> problem(eta);
  
 //right hand side
  auto flambda = [&](const auto& e, const auto& x)
   {return problem.f(e,x);};
  auto f = Dune::PDELab::
    makeGridFunctionFromCallable(gv,flambda);
 
  //Dirichlet extension
  auto glambda = [&](const auto& e, const auto& x)
    {return problem.g(e,x);};
  auto g = Dune::PDELab::
    makeGridFunctionFromCallable(gv,glambda);
  
  //Dirichlet type of BC
  auto blambda = [&](const auto& i, const auto& x)
    {return problem.b(i,x);};
  auto b = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,blambda);
  
// Make grid function space
  
  typedef Dune::PDELab::ConformingDirichletConstraints CON; 
  typedef Dune::PDELab::istl::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("Pk");

  // Assemble constraints
  typedef typename GFS::template
    ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;

  // A coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value

  // Make a grid function out of it
  typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  ZDGF zdgf(gfs,z);

  // Fill the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,z);

// make vector consistent NEW IN PARALLEL //DIFF//////////////////////////////////////////////////
  Dune::PDELab::ISTL::ParallelHelper<GFS> helper(gfs);
  helper.maskForeignDOFs(z);
  Dune::PDELab::AddDataHandle<GFS,Z> adddh(gfs,z);
  if (gfs.gridView().comm().size()>1)
    gfs.gridView().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
  //////////////////////////////////////////////////////////////////////////////////////////////////

 // Make a local operator
  typedef PoissonPk< decltype(f),FEM > LOP;
  LOP lop(f,fem.find(*gv.template begin<0>()));

 /* Make a local operator
  typedef PoissonPk<problem<RF>,FEM> LOP;
  LOP lop(problem);*/

  // Make a global operator
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  //MBE mbe(1<<(dim+1)); // guess nonzeros per row
  int degree = ptree.get("fem.degree",(int)1);
  MBE mbe((int)pow(1+2*degree,dim));
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,  /* ansatz and test space */
    LOP,      /* local operator */
    MBE,      /* matrix backend */
    RF,RF,RF, /* domain, range, jacobian field type*/
    CC,CC     /* constraints for ansatz and test space */
    > GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // Select a linear solver backend NEW IN PARALLEL//////////DIFF
  typedef Dune::PDELab::ISTLBackend_CG_AMG_SSOR<GO> LS;
  int verbose=0;//DIFF
  if (gfs.gridView().comm().rank()==0) verbose=1;//DIFF
  LS ls(gfs,100,verbose);//DIFF


  // Select a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS ;
  //LS ls(5000,true);
  

 // typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
  //LS ls(100,3);

  // Assemble and solve linear problem
  typedef Dune::PDELab::
    StationaryLinearProblemSolver<GO,LS,Z> SLP;
  SLP slp(go,ls,z,1e-10,1);
  slp.apply(); // here all the work is done! 

 


  //Z w(gfs); // Lagrange interpolation of exact solution
 // Dune::PDELab::interpolate(g,gfs,w);
  // ZDGF wdgf(gfs,w);
  
 
  // Write VTK output file

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming); 
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  //vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
    //                                     VTKF(wdgf,"exact")));
  vtkwriter.write(ptree.get("output.filename","output"),
                  Dune::VTK::appendedraw);


 //Here, the calculation of l2 norm error starts

// auto l2norm = l2difference(wdgf,zdgf,1);
 //std::cout  << "error is =" << l2norm << std::endl;
 
  
        
 
}
