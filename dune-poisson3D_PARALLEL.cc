// -*- tab-width: 4; indent-tabs-mode: nil -*-
// Beware of line number changes, they may corrupt docu!
/** \file
    \brief Solve Poisson equation with P1 conforming finite elements
*/
// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh> //DIFF
#include <dune/grid/common/grid.hh>//added
#include<dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/parmetisgridpartitioner.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// include all components making up this tutorial
//#include"poissonp1.hh"
#include"PoissonPk.hh"
#include"driver.hh"
//#include "problemCube.hh"
//includes for l2 norm calculation in the driver.hh
#include<cmath>
#include<vector>
#include <parmetis.h>
//#undef HAVE_DUNE_ALUGRID
//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    // Maybe initialize Mpi
    Dune::MPIHelper&
      mpiHelper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout << "Parallel code run on "
                << mpiHelper.size() << " process(es)" << std::endl;

    // open ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("dune-poisson3D_PARALLEL.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // read parameteres from ini file
    const int dim = ptree.get<int>("grid.dim");
    const int refinement = ptree.get<int>("grid.refinement");
    std::string gridmanager = ptree.get("grid.manager","yasp");//added
    const int degree = ptree.get<int>("fem.degree");
    //get partition type 
     std::string type = ptree.get("partition.type","dune");


    ////////////////ADDED//////////////////////////////////


// use UG grid if available and selected
 
   if (dim==2 && gridmanager=="ug")
      {
#if  HAVE_UG
  typedef Dune::UGGrid<2> Grid;

  std::string filename = ptree.get("grid.twod.filename","unitsquare.msh");
 
  Dune::GridFactory<Grid> factory;
  Dune::GmshReader<Grid>::read(factory,filename,true,true);
  std::shared_ptr<Grid> gridp(factory.createGrid());
  
  Dune::Timer timer;
  gridp->globalRefine(refinement);
  std::cout << "Time for mesh refinement " << timer.elapsed()
                  << " seconds" << std::endl;
  if(type=="dune"){
  gridp->loadBalance();}
  
  typedef Grid::LeafGridView GV;
  typedef Grid::ctype DF;
  GV gv=gridp->leafGridView();

   /*if(type=="parMetis"){

  // Create initial partitioning using ParMETIS
   std::vector<unsigned int> part
    = ParMetisGridPartitioner<GV>::partition(gv, mpiHelper);

    gridp->loadBalance(part,0);}*/


  if (degree==1) {
         typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
   }
   if (degree==2) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
   }
   if (degree==3) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,3> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
   }

#else
        std::cout << "You selected ug as grid manager but ug was not found during installation" << std::endl;
#endif
     }
    if (dim==3 && gridmanager=="ug")
      {
#if HAVE_UG
        typedef Dune::UGGrid<3> Grid;
        std::string filename = ptree.get("grid.threed.filename","nut.msh");
        

  Dune::GridFactory<Grid> factory;
  Dune::GmshReader<Grid>::read(factory,filename,true,true);
  std::shared_ptr<Grid> gridp(factory.createGrid());
  
  Dune::Timer timer;
  gridp->globalRefine(refinement);
  std::cout << "Time for mesh refinement " << timer.elapsed()
                  << " seconds" << std::endl;
   if(type=="dune"){
     gridp->loadBalance();}
  
  typedef Grid::LeafGridView GV;
  typedef Grid::ctype DF;
  GV gv=gridp->leafGridView();

/*  if(type=="parMetis"){

  // Create initial partitioning using ParMETIS
   std::vector<unsigned int> part
    = ParMetisGridPartitioner<GV>::partition(gv, mpiHelper);

    gridp->loadBalance(part,0);}*/

        if (degree==1) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==2) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==3) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,3> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
#else
        std::cout << "You selected ug as grid manager but ug was not found during installation" << std::endl;
#endif
      }

if (dim==2 && gridmanager=="alu")
      {
#if HAVE_DUNE_ALUGRID
        typedef Dune::ALUGrid<2,2,Dune::simplex,
                              Dune::nonconforming> Grid;
        std::string filename = ptree.get("grid.twod.filename","unitsquare.msh");

  Dune::GridFactory<Grid> factory;
  Dune::GmshReader<Grid>::read(factory,filename,true,true);
  std::shared_ptr<Grid> gridp(factory.createGrid());
  
  Dune::Timer timer;
  gridp->globalRefine(refinement);
  std::cout << "Time for mesh refinement " << timer.elapsed()
                  << " seconds" << std::endl;
  if(type=="dune"){
   gridp->loadBalance();}
  
  typedef Grid::LeafGridView GV;
  typedef Grid::ctype DF;
  GV gv=gridp->leafGridView();

 /* if(type=="parMetis"){

  // Create initial partitioning using ParMETIS
   std::vector<unsigned int> part
    = ParMetisGridPartitioner<GV>::partition(gv, mpiHelper);

    gridp->loadBalance(part,0);}*/

        if (degree==1) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==2) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==3) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,3> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
#else
        std::cout << "You selected alugrid as grid manager but alugrid was not found during installation" << std::endl;
#endif
      }
    if (dim==3 && gridmanager=="alu")
      {
#if HAVE_DUNE_ALUGRID
        typedef Dune::ALUGrid<3,3,Dune::simplex,
                              Dune::nonconforming> Grid;
        std::string filename = ptree.get("grid.threed.filename","nut.msh");
        

  Dune::GridFactory<Grid> factory;
  Dune::GmshReader<Grid>::read(factory,filename,true,true);
  std::shared_ptr<Grid> gridp(factory.createGrid());
  
  Dune::Timer timer;
  gridp->globalRefine(refinement);
  std::cout << "Time for mesh refinement " << timer.elapsed()
                  << " seconds" << std::endl;
  if(type=="dune"){
  gridp->loadBalance();}
  
  typedef Grid::LeafGridView GV;
  typedef Grid::ctype DF;
  GV gv=gridp->leafGridView();

/*  if(type=="parMetis"){

  // Create initial partitioning using ParMETIS
   std::vector<unsigned int> part
    = ParMetisGridPartitioner<GV>::partition(gv, mpiHelper);

    gridp->loadBalance(part,0);}*/

        if (degree==1) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==2) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==3) {
          typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,3> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
#else
        std::cout << "You selected alugrid as grid manager but alugrid was not found during installation" << std::endl;
#endif
        }

    ///////////////ADDED///////////////////////////////////
     
   // YaspGrid section
    if (dim==2 && gridmanager=="yasp" )
      {
        const int dim=2;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<DF,dim> L;
        L[0] = ptree.get("grid.structured.LX",(double)1.0);
        L[1] = ptree.get("grid.structured.LY",(double)1.0);
        std::array<int,dim> N; //grid resolution
        N[0] = ptree.get("grid.structured.NX",(int)25);
        N[1] = ptree.get("grid.structured.NY",(int)25);
        std::bitset<dim> B(false);//diff //no periodic in either  dimensions
        int overlap=1;//diff  //thickness of the overlap layer
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new     Grid(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication())); //DIFF
        gridp->refineOptions(false); // keep overlap in cells   //DIFF
        gridp->globalRefine(refinement);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();

        if (degree==1) {
          typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }

        if (degree==2) {
          typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
      }
    if (dim==3 && gridmanager=="yasp")
      {
        const int dim=3;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<DF,dim> L;
        L[0] = ptree.get("grid.structured.LX",(double)1.0);
        L[1] = ptree.get("grid.structured.LY",(double)1.0);
        L[2] = ptree.get("grid.structured.LZ",(double)1.0);
        std::array<int,dim> N;
        N[0] = ptree.get("grid.structured.NX",(int)25);
        N[1] = ptree.get("grid.structured.NY",(int)25);
        N[2] = ptree.get("grid.structured.NZ",(int)25);
        std::bitset<dim> B(false);//diff
        int overlap=1;//diff
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N,B,overlap,Dune::MPIHelper::getCollectiveCommunication())); //DIFF
        gridp->refineOptions(false); // keep overlap in cells       //DIFF
        gridp->globalRefine(refinement);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
        if (degree==1) {
          typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,1> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
        if (degree==2) {
          typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,double,2> FEM;
          FEM fem(gv);
          driver(gv,fem,ptree);
        }
      }
  }


catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
