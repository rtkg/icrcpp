#include "../include/search_zones.h"
#include "../include/debug.h"
#include "../include/grasp.h"
#include "assert.h"
#include <libqhullcpp/QhullVertexSet.h>


#include <iomanip>
#include <sys/time.h>
#include <time.h>

namespace ICR
{
#define KAPPA 100000
  CasADi::FX generateCodeAndCompile(CasADi::FX fcn, const std::string& name)
  {
    std::cout << "Generating code for " << name << std::endl;

    // // Convert to an SXFunction (may or may not improve efficiency)
    // if( is_a<MXFunction>(fcn)){
    //   fcn = SXFunction(shared_cast<MXFunction>(fcn));
    //   fcn.init();
    // }

    // Generate C code
    fcn.generateCode(name + ".c");

    // Compilation command
    std::string compile_command = "gcc -fPIC -shared -pipe " + name + ".c -o " + name + ".so";

    // Compile the c-code
    int flag = system(compile_command.c_str());
    casadi_assert_message(flag==0, "Compilation failed");

    // Load the generated function for evaluation
    CasADi::ExternalFunction fcn_e("./" + name + ".so");
    return fcn_e;
  }
  //-------------------------------------------------------------------
  PrimitiveSearchZone::PrimitiveSearchZone(){hyperplane_ids_.clear();}
  //-------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream, PrimitiveSearchZone const& psz)
  {
    stream <<'\n'<<"PRIMITIVE SEARCH ZONE: "<<'\n'
	   <<"Hyperplane id's: ";
       
    for(uint i=0;i<psz.hyperplane_ids_.size();i++)
      stream<<psz.hyperplane_ids_[i]<<" ";

    stream<<'\n'<<'\n';

    return stream;
  }
  //-------------------------------------------------------------------
  SearchZones::SearchZones() : num_search_zones_(0),search_zones_computed_(false)
  {
    search_zones_.clear();
  }
  //-------------------------------------------------------------------
  SearchZones::SearchZones(const GraspPtr grasp) : grasp_(grasp), num_search_zones_(0),search_zones_computed_(false)
  {
    assert(grasp_->isInitialized());
    assert(grasp_->getGWS()->containsOrigin());
    search_zones_.clear();
  }
  //-------------------------------------------------------------------
  SearchZones::SearchZones(SearchZones const& src) :  grasp_(src.grasp_), tws_(src.tws_), search_zones_(src.search_zones_),
						      num_search_zones_(src.num_search_zones_), search_zones_computed_(src.search_zones_computed_),map_vertex2finger_(src.map_vertex2finger_),
                                                      hyperplane_normals_(src.hyperplane_normals_),hyperplane_offsets_(src.hyperplane_offsets_){}
  //-------------------------------------------------------------------
  SearchZones& SearchZones::operator=(SearchZones const& src)
  {
    if (this !=&src)
      {
	grasp_=src.grasp_;
	tws_=src.tws_;
	search_zones_=src.search_zones_;
	num_search_zones_=src.num_search_zones_;
	search_zones_computed_=src.search_zones_computed_;
	map_vertex2finger_=src.map_vertex2finger_;
	hyperplane_normals_=src.hyperplane_normals_;
	hyperplane_offsets_=src.hyperplane_offsets_;
      }
    return *this;
  }
  //-------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& stream,SearchZones const& sz)
  {
    stream <<'\n'<<"SEARCH ZONES: "<<'\n'
	   <<"Number of patch search zones: " << sz.num_search_zones_<<'\n'
	   <<"Search zones computed: "<<std::boolalpha<<sz.search_zones_computed_<<'\n'
	   <<"Vertex-finger map: "<<sz.map_vertex2finger_<<'\n';

    return stream;
  }
  //-------------------------------------------------------------------
  SearchZones::~SearchZones(){clear();}
  //-------------------------------------------------------------------
  SearchZone const* SearchZones::getSearchZone(uint finger_id)const{return search_zones_.at(finger_id);}
  //-------------------------------------------------------------------
  uint SearchZones::getNumSearchZones()const{return num_search_zones_;}
  //-------------------------------------------------------------------
  bool SearchZones::searchZonesComputed()const{return search_zones_computed_;}
  //-------------------------------------------------------------------
  void SearchZones::computeShiftedHyperplanes()
  {
    uint K=grasp_->getGWS()->getDimension();
    uint H=grasp_->getGWS()->num_facets_;
    hyperplane_normals_.resize(H,K);
    hyperplane_offsets_.resize(H);
    facetT* curr_f=grasp_->getGWS()->conv_hull_.beginFacet().getFacetT();
 
    if(tws_->getWrenchSpaceType() == Spherical)
      {
	double offset=tws_->getOcInsphereRadius();
	hyperplane_offsets_.fill(offset); 
     
	for(uint h=0; h < H;h++)
	  {
	    assert(offset <= (-curr_f->offset)); //Make sure that GWSinit contains the TWS and the hyperplanes are shifted inwards
	    //The direction of the normal is switched to inward-pointing in order to be consistent with the positive offset

	    //  hyperplane_normals_.row(h)=(-Eigen::Map<Eigen::Matrix<double,1,6> >(curr_f->normal));
	    hyperplane_normals_.row(h)=-Eigen::Map<Eigen::RowVectorXd>(curr_f->normal,K);
	    curr_f=curr_f->next;
	  }
      }

    else if(tws_->getWrenchSpaceType() == Discrete)
      {
	//Should assert that the TWS contains the origin ...
	Eigen::MatrixXd TW(K,dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get())->getNumWrenches());
        doubleArrayToEigenMatrix(dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get())->getWrenches().get(),TW);

      	for(uint h=0; h < H;h++)
	  {
	    Eigen::RowVectorXd normal=-Eigen::Map<Eigen::RowVectorXd>(curr_f->normal,K);
            double e = (-curr_f->offset);
	    assert(e > 0); //just in case ...

	    double e_shift=-(normal*TW).minCoeff();
	    //Should provide some tolerances for the checks below ...
	    if (e_shift < 0)
	      {
		std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - the TWS does not contain the origin! Exiting ..."<<std::endl;
		exit(0);
             }
	    if (e_shift > e)
	      {
		std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - the prototype grasp's GWS does not contain the TWS! Exiting ..."<<std::endl;
		exit(0);
              }


	    hyperplane_offsets_(h)=e_shift;
	    hyperplane_normals_.row(h)=normal;
	    // std::cout<<hyperplane_normals_.row(h)<<" "<<hyperplane_offsets_(h)<<std::endl;
           curr_f=curr_f->next;
          }

      }
    else
      std::cout<<"Error in SearchZones::computeShiftedHyperplanes() - invalid Task Wrench Space type!"<<std::endl;
  }
  //---------------------------------------------------------------------
  void SearchZones::setTWSNLPHyperplanes(const CasADi::IpoptSolver& nlp_solver,const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list)
  {
    uint H=Wh_i_list.size();
    uint L=Wh_i_list[0][0]->cols(); //dimension
    uint K=Wh_i_list[0][0]->rows(); //dimension
    CasADi::DMatrix solution=nlp_solver.output("x");

    hyperplane_normals_.resize(H,L);
    hyperplane_offsets_.resize(H);
      
    Eigen::RowVectorXd normal(1,K);
    uint h_start_x_ind=0;
    for (uint h=0; h<H;h++)
      {
	uint Ih=Wh_i_list[h].size();
	uint nVh=K+1+Ih*(L+1);

	for(uint j=0; j<K; j++)
	  normal(j)=solution(h_start_x_ind+j).toScalar();

	double de=1/normal.norm();
	normal=normal*de;

	double offset=solution(h_start_x_ind+K).toScalar()*de;
	assert(offset >= 0);
	hyperplane_normals_.row(h)=normal;
	hyperplane_offsets_(h)=offset;


	h_start_x_ind+=nVh;
      }

  }
  //-------------------------------------------------------------------------------------------------
  void SearchZones::initialSolutionTWSNLP(CasADi::DMatrix& initial_solution, const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list)
  {
    //wrenches in Ph_list and Wh_i_list are stored as column - matrices

    uint H=Wh_i_list.size();
    uint nV=0; //number of decision variables
    uint nC=0; //number of constraints
    uint K=Wh_i_list[0][0]->rows(); //dimension
    uint L=Wh_i_list[0][0]->cols(); //number of wrenches in a wrench cone

    DiscreteTaskWrenchSpace* tws=dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get());
    uint T=tws->getNumWrenches();

    for (uint h=0; h<H;h++)
      {
	nV+=Wh_i_list[h].size()*(L+1)+2;
	nC+=2*Wh_i_list[h].size()+T;
      }


    //Extract the task wrenches and put them in a CasADi::DMatrix
    CasADi::DMatrix TW(K,T);
    uint count=0;
    for (int j=0; j<T; j++)
      for (int i=0; i<K; i++)
	{
	  TW(i,j)=tws->getWrenches().get()[count];
	  count++;
	}

    //decision variable vector
    CasADi::SXMatrix x = CasADi::ssym("x",nV); 
    CasADi::DMatrix x0=CasADi::DMatrix::zeros(nV,1);

    //Objective
    CasADi::SXMatrix f(1,1); 

    //Constraints
    CasADi::SXMatrix g=CasADi::SXMatrix::zeros(nC,1); 
    CasADi::DMatrix ub_g=CasADi::DMatrix::inf(nC,1);
    CasADi::DMatrix lb_g=CasADi::DMatrix::inf(nC,1)*(-1);

    //Bounds
    CasADi::DMatrix ub_x=CasADi::DMatrix::inf(nV,1);
    CasADi::DMatrix lb_x=CasADi::DMatrix::inf(nV,1)*(-1);

    //Index Sets for the decision variables
    double eh_ind;
    double deh_ind;
    std::vector<int> lmbdh_ind;
    std::vector<std::vector<int> > lmbdh_i_ind;
    std::vector<int> zh_i_ind;

    //Index Sets for the constraints
    std::vector<int> tc_ind(T);
    std::vector<int> eqc_ind;
    std::vector<int> sc_ind;

    uint h_start_x_ind=0;
    uint h_start_c_ind=0;

    facetT* curr_f=grasp_->getGWS()->conv_hull_.beginFacet().getFacetT();
    for (uint h=0; h<H;h++)
      {
	//================FORM INDEX SETS========================================
	uint Ih=Wh_i_list[h].size();
	uint nVh=2+Ih*(L+1);

	//Build the h-th index lists for nh, eh, lbmd_h_i and zh_i for eache Hyperplane h
	lmbdh_ind.resize(Ih*L);
	lmbdh_i_ind.resize(Ih);
	zh_i_ind.resize(Ih);

   
	eh_ind=h_start_x_ind;
	deh_ind=h_start_x_ind+1;

	for (int i=0; i<Ih; i++)
	  {
	    lmbdh_i_ind[i].resize(L);
	    for(int l=0; l<L; l++)
	      lmbdh_i_ind[i][l]=L*i+l+h_start_x_ind+2;

	    zh_i_ind[i]=h_start_x_ind+2+L*Ih+i;
	  }

	for (uint i=0; i<Ih*L; i++)
	  lmbdh_ind[i]=h_start_x_ind+2+i;

	//Build the h-th index lists for the constriants for each Hyperplane h
	uint nCh=2*Ih+T;
	eqc_ind.resize(Ih);
	sc_ind.resize(Ih);

	for (uint t=0; t<T;t++)
	  tc_ind[t]=h_start_c_ind+t;

	for (uint i=0; i<Ih;i++)
	  {
	    eqc_ind[i]=h_start_c_ind+T+i;
	    sc_ind[i]=h_start_c_ind+T+Ih+i;
	  }

	h_start_x_ind+=nVh;
	h_start_c_ind+=nCh;

	//================OBJECTIVE========================================
	f=f-x(deh_ind)-KAPPA*mul(x(zh_i_ind).trans(),CasADi::SXMatrix::ones(Ih,1));

	//==================== BOUNDS ========================================
	//de's have to be positive
	lb_x(deh_ind)=CasADi::DMatrix::zeros(1,1);

	//Slack variables have to be negative
	ub_x(zh_i_ind)=CasADi::DMatrix::zeros(Ih,1);

	//Lambdas have to be  positive
	for (int i=0; i<Ih;i++)
	  lb_x(lmbdh_i_ind[i])=CasADi::DMatrix::zeros(L,1);

	//====================HARD CONSTRAINTS ========================================

	//Constraints involving the task wrenches
	//get the h-th normal of the prototype grasp's GWS
	CasADi::DMatrix ah(K,1); 
	for (int k=0; k<K;k++)
          ah(k)=curr_f->normal[k]*(-1);

	g(tc_ind) = mul(x(eh_ind),-CasADi::SXMatrix::ones(T,1));
	ub_g(tc_ind) = mul(TW.trans(),ah);

	//Lambdas have to sum up to one
	//Build the Equality matrix
	std::vector<int> col_ind(L);
	CasADi::SXMatrix E=CasADi::SXMatrix::zeros(Ih,Ih*L);
	for (uint i=0; i<Ih;i++)
	  {
	    for (uint j=0; j<L;j++)
	      col_ind[j]=i*L+j;
	     
	    E(i,col_ind)=CasADi::SXMatrix::ones(1,L);
	  }
	//set the Equality constraints
	g(eqc_ind)=mul(E,x(lmbdh_ind));

	ub_g(eqc_ind)=CasADi::DMatrix::ones(Ih,1);
	lb_g(eqc_ind)=CasADi::DMatrix::ones(Ih,1);

	//==================== SOFT CONSTRAINTS ======================================
	for (uint i=0; i<Ih; i++)
	  {
	    //Extract the i-th Wrench cone and store it in a CasADi::DMatrix
	    CasADi::SXMatrix Whi(Wh_i_list[h][i]->rows(),Wh_i_list[h][i]->cols());
	    for (uint r=0; r<Wh_i_list[h][i]->rows(); r++)
	      for (uint c=0;c<Wh_i_list[h][i]->cols() ; c++)
		Whi(r,c)=(*Wh_i_list[h][i])(r,c);


	    g(sc_ind[i])=mul(CasADi::SXMatrix(ah.trans()), mul(Whi,x(lmbdh_i_ind[i]))) +x(eh_ind)+x(deh_ind)+x(zh_i_ind[i]);
	    ub_g(sc_ind[i])=CasADi::DMatrix::zeros(1,1);
	  }

	curr_f=curr_f->next;
      }		  

    //  std::cout<<"f: "<<std::endl<<f<<std::endl;
    //  std::cout<<"g: "<<std::endl<<g<<std::endl;
    //  std::cout<<"ub_g: "<<std::endl<<ub_g<<std::endl;
    //  std::cout<<"lb_g: "<<std::endl<<lb_g<<std::endl;
    //  std::cout<<"ub_x: "<<std::endl<<ub_x<<std::endl;
    //  std::cout<<"lb_x: "<<std::endl<<lb_x<<std::endl;

    CasADi::SXFunction lp(nlpIn("x",x),nlpOut("f",f,"g",g));
    CasADi::IpoptSolver lp_solver(lp);
    //  CasADi::QPSolver qp_solver();


    lp_solver.setOption("print_level",0);
    lp_solver.setOption("linear_solver","ma97");
    lp_solver.setOption("mu_strategy","adaptive");
    lp_solver.init();

    lp_solver.setInput( x0, "x0");
    lp_solver.setInput(lb_x, "lbx");
    lp_solver.setInput(ub_x, "ubx");
    lp_solver.setInput(lb_g, "lbg");
    lp_solver.setInput(ub_g, "ubg");

    lp_solver.evaluate();

    // std::cout << std::setw(30) << "Initial Objective: " << lp_solver.output("f").getDescription() << std::endl;
    // std::cout << std::setw(30) << "Initial Primal solution: " << lp_solver.output("x").getDescription() << std::endl;


    //Extract the initial vectors from the solution
    initial_solution=CasADi::DMatrix();
    curr_f=grasp_->getGWS()->conv_hull_.beginFacet().getFacetT();
    h_start_x_ind=0;
    CasADi::DMatrix solution=lp_solver.output("x");
    for (uint h=0; h<H; h++)
      {
	uint Ih=Wh_i_list[h].size();
	uint nVh=2+Ih*(L+1);
    
	eh_ind=h_start_x_ind;
	deh_ind=h_start_x_ind+1;
	std::vector<int> lmbdh_zh_ind(Ih*(L+1));
	for (int j=0; j<lmbdh_zh_ind.size(); j++)
	  lmbdh_zh_ind[j]=h_start_x_ind+2+j;
  
	//get the h-th normal
	CasADi::DMatrix ah=CasADi::DMatrix::zeros(K,1);
	for (int k=0; k<K;k++)
          ah(k)=curr_f->normal[k]*(-1);


	//HAAAAACKKKKK to avoid division by zero ... should rethink that ...
	double deh=solution(deh_ind).toScalar();
	if(deh > 1e-5)
	  {
	    initial_solution=CasADi::vertcat(initial_solution, ah/deh); 
	    initial_solution=CasADi::vertcat(initial_solution, solution(eh_ind)/deh); 
	  }
	else
	  {
	    //	std::cout<<"Violation: deh="<<deh<<" eh="<<solution(eh_ind)<<std::endl;
	    initial_solution=CasADi::vertcat(initial_solution, ah*(0)); 
	    initial_solution=CasADi::vertcat(initial_solution, solution(eh_ind)*0); 
	  }
	initial_solution=CasADi::vertcat(initial_solution, solution(lmbdh_zh_ind));

 
	h_start_x_ind+=nVh;
	curr_f=curr_f->next;
      }
 
  }
  //-------------------------------------------------------------------
  void SearchZones::createDiscreteTWSNLPSolver(const std::vector<Eigen::MatrixXd*>& Ph_list,const std::vector<std::vector<Eigen::MatrixXd* > >& Wh_i_list,CasADi::IpoptSolver& nlp_solver)
  {
    //wrenches in Ph_list and Wh_i_list are stored as column - matrices

    uint H=Ph_list.size();
    uint nV=0; //number of decision variables
    uint nC=0; //number of constraints
    uint K=Wh_i_list[0][0]->rows(); //dimension
    uint L=Wh_i_list[0][0]->cols(); //number of wrenches in a wrench cone

    DiscreteTaskWrenchSpace* tws=dynamic_cast<DiscreteTaskWrenchSpace*>(tws_.get());
    uint T=tws->getNumWrenches();

    for (uint h=0; h<H;h++)
      {
	nV+=Wh_i_list[h].size()*(L+1)+K+1;
	nC+=2*Wh_i_list[h].size()+T+Ph_list[h]->cols();
      }


    //Extract the task wrenches and put them in a CasADi::SXMatrix
    CasADi::SXMatrix TW(K,T);
    uint count=0;
    for (int j=0; j<T; j++)
      for (int i=0; i<K; i++)
	{
	  TW(i,j)=tws->getWrenches().get()[count];
	  count++;
	}


    // Utility for timing selected parts of the code - uncomment below and put the code to be timed at the marked location
    struct timeval start, end;
    double c_time;
    // gettimeofday(&start,0);
    //
    //  Put code to be timed here...
    //  
    // gettimeofday(&end,0);
    // c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
    // std::cout<<"Computation time: "<<c_time<<" s"<<std::endl;

    //decision variable vector & initial solution
    CasADi::SXMatrix x = CasADi::ssym("x",nV); 
    CasADi::DMatrix x0=CasADi::DMatrix::zeros(nV,1); 
    gettimeofday(&start,0);
      initialSolutionTWSNLP(x0,Wh_i_list);
    gettimeofday(&end,0);
    c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
    std::cout<<"Computation time comute initial solution: "<<c_time<<" s"<<std::endl;

    //Objective
    CasADi::SXMatrix f(1,1); 

    //Constraints
    CasADi::SXMatrix g=CasADi::SXMatrix::zeros(nC,1); 
    CasADi::DMatrix ub_g=CasADi::DMatrix::inf(nC,1);
    CasADi::DMatrix lb_g=CasADi::DMatrix::inf(nC,1)*(-1);

    //Bounds
    CasADi::DMatrix ub_x=CasADi::DMatrix::inf(nV,1);
    CasADi::DMatrix lb_x=CasADi::DMatrix::inf(nV,1)*(-1);

    //Index Sets for the decision variables
    std::vector<int> nh_ind(K);
    double eh_ind;
    std::vector<std::vector<int> > lmbdh_i_ind;
    std::vector<int> lmbdh_ind;
    std::vector<int> zh_i_ind;

    //Index Sets for the constraints
    std::vector<int> phc_ind;
    std::vector<int> tc_ind(T);
    std::vector<int> eqc_ind;
    std::vector<int> sc_ind;

    uint h_start_x_ind=0;
    uint h_start_c_ind=0;


    //std::cout<<"Number of hyperplanes to fit: "<<H<<std::endl;

    gettimeofday(&start,0);
    for (uint h=0; h<H;h++)
      {

	//===============FORM INDEX SETS========================================
	uint Ih=Wh_i_list[h].size();
	uint n_Ph=Ph_list[h]->cols();


	//Build the h-th index lists for nh, eh, lbmd_h_i and zh_i for each Hyperplane h
	uint nVh=K+1+Ih*(L+1);
	lmbdh_i_ind.resize(Ih);
	lmbdh_ind.resize(Ih*L);
	zh_i_ind.resize(Ih);

	for (uint j=0; j<K; j++)
	  nh_ind[j]=h_start_x_ind+j;

	eh_ind=h_start_x_ind+K;

	for (uint i=0; i<Ih; i++)
	  {
	    lmbdh_i_ind[i].resize(L);
	    for(uint l=0; l<L; l++)
	      lmbdh_i_ind[i][l]=L*i+l+h_start_x_ind+K+1;

	    zh_i_ind[i]=h_start_x_ind+K+1+L*Ih+i;
	  }

	for (uint i=0; i<Ih*L; i++)
	  lmbdh_ind[i]=h_start_x_ind+K+1+i;

	//Build the h-th index lists for the constriants for each Hyperplane h
	uint nCh=2*Ih+T+n_Ph;  
	phc_ind.resize(n_Ph);
	eqc_ind.resize(Ih);
	sc_ind.resize(Ih);


	for (uint p=0; p<n_Ph;p++)
	  phc_ind[p]=h_start_c_ind+p;

	for (uint t=0; t<T;t++)
	  tc_ind[t]=h_start_c_ind+n_Ph+t;

	for (uint i=0; i<Ih;i++)
	  {
	    eqc_ind[i]=h_start_c_ind+n_Ph+T+i;
	    sc_ind[i]=h_start_c_ind+n_Ph+T+Ih+i;
	  }

	h_start_x_ind+=nVh;
	h_start_c_ind+=nCh;

	//================OBJECTIVE========================================
	f=f+0.5*mul(x(nh_ind).trans(),x(nh_ind))-KAPPA*mul(x(zh_i_ind).trans(),CasADi::SXMatrix::ones(Ih,1));

	//================HARD CONSTRAINTS========================================
	//Extract the wrenches of the prototype grasp and store them in a CasADi::SXMatrix
	CasADi::SXMatrix Ph(Ph_list[h]->rows(),Ph_list[h]->cols());
	for (uint i=0; i<Ph_list[h]->rows(); i++)
	  for (uint j=0; j<Ph_list[h]->cols(); j++)
	    Ph(i,j)=(*Ph_list[h])(i,j);

	//Constraints involving the prototype grasps's wrenches Ph_i
	g(phc_ind)=mul(Ph.trans(),x(nh_ind))+mul(x(eh_ind),CasADi::SXMatrix::ones(K,1));
	ub_g(phc_ind)=CasADi::DMatrix(n_Ph,1,-1);

	//Constraints involving the task wrenches
	g(tc_ind) = -mul(TW.trans(),x(nh_ind))-mul(x(eh_ind),CasADi::SXMatrix::ones(T,1));
	ub_g(tc_ind)=CasADi::DMatrix(T,1,0); //0 for close fit, -1 for middle fit

	//Lambdas have to sum up to one
	//Build the Equality matrix
	std::vector<int> col_ind(L);
	CasADi::SXMatrix E=CasADi::SXMatrix::zeros(Ih,Ih*L);
	for (uint i=0; i<Ih;i++)
	  {
	    for (uint j=0; j<L;j++)
	      col_ind[j]=i*L+j;
	     
	    E(i,col_ind)=CasADi::SXMatrix::ones(1,L);
	  }

	g(eqc_ind)=mul(E,x(lmbdh_ind));
	ub_g(eqc_ind)=CasADi::DMatrix::ones(Ih,1);
	lb_g(eqc_ind)=CasADi::DMatrix::ones(Ih,1);


	//=================SOFT CONSTRAINTS========================================
	for (uint i=0; i<Ih; i++)
	  {
	    //Extract the i-th Wrench cone and store it in a CasADi::SXMatrix
	    CasADi::SXMatrix Whi(Wh_i_list[h][i]->rows(),Wh_i_list[h][i]->cols());
	    for (uint r=0; r<Wh_i_list[h][i]->rows(); r++)
	      for (uint c=0;c<Wh_i_list[h][i]->cols() ; c++)
		Whi(r,c)=(*Wh_i_list[h][i])(r,c);

	    g(sc_ind[i])=mul(x(nh_ind).trans(),mul(Whi,x(lmbdh_i_ind[i])))+x(eh_ind)+x(zh_i_ind[i]);
	    ub_g(sc_ind[i])=-CasADi::DMatrix::ones(1,1);

	  }
 
	//==================== BOUNDS ========================================
	//Slack variables have to be negative
	ub_x(zh_i_ind)=CasADi::DMatrix::zeros(Ih,1);

	//Lambdas have to be  positive
	for (uint i=0; i<Ih;i++)
	  lb_x(lmbdh_i_ind[i])=CasADi::DMatrix::zeros(L,1);

      }

    //  std::cout<<"f: "<<std::endl<<f<<std::endl;
    //  std::cout<<"g: "<<std::endl<<g<<std::endl;
    // std::cout<<"ub_g: "<<std::endl<<ub_g<<std::endl;
    // std::cout<<"lb_g: "<<std::endl<<lb_g<<std::endl;
    // std::cout<<"ub_x: "<<std::endl<<ub_x<<std::endl;
    // std::cout<<"lb_x: "<<std::endl<<lb_x<<std::endl;

    gettimeofday(&end,0);
    c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
    std::cout<<"Computation time Formulate Problem: "<<c_time<<" s"<<std::endl<<std::endl;

    gettimeofday(&start,0);

    bool codegen=false;
    if(codegen)
      {
	//////////////////CODEGEN//////////////////////////////////////////
	// NLP function
	CasADi::FX nlp = CasADi::SXFunction(CasADi::nlpIn("x",x),CasADi::nlpOut("f",f,"g",g));
	nlp.init();
	// Gradient of the objective
	CasADi::FX grad_f = nlp.gradient("x","f");
	grad_f.init();
	// Jacobian of the constraints
	CasADi::FX jac_g = nlp.jacobian("x","g");
	jac_g.init();
	// Hessian of the lagrangian
	CasADi::FX grad_lag = nlp.derivative(0,1);
	CasADi::FX hess_lag = grad_lag.jacobian(CasADi::NL_X,CasADi::NL_NUM_OUT+CasADi::NL_X,false,true);
	hess_lag.init();
	// Codegen and compile
	nlp = generateCodeAndCompile(nlp,"nlp");
	grad_f = generateCodeAndCompile(grad_f,"grad_f");
	jac_g = generateCodeAndCompile(jac_g,"jac_g");
	hess_lag = generateCodeAndCompile(hess_lag,"hess_lag");
	// Create an NLP solver passing derivative information
	nlp_solver=CasADi::IpoptSolver(nlp);
	nlp_solver.setOption("grad_f",grad_f);
	nlp_solver.setOption("jac_g",jac_g);
	nlp_solver.setOption("hess_lag",hess_lag);
      }
    else
      {
	//////////////////NO CODEGEN//////////////////////////////////////////
	CasADi::SXFunction nlp(CasADi::nlpIn("x",x),CasADi::nlpOut("f",f,"g",g));
	nlp_solver=CasADi::IpoptSolver(nlp);
      }

    nlp_solver.setOption("print_level",0);
    nlp_solver.setOption("linear_solver","ma97");
    nlp_solver.setOption("mu_strategy","adaptive");
    nlp_solver.init();

    nlp_solver.setInput( x0, "x0");
    nlp_solver.setInput(lb_x, "lbx");
    nlp_solver.setInput(ub_x, "ubx");
    nlp_solver.setInput(lb_g, "lbg");
    nlp_solver.setInput(ub_g, "ubg");

    gettimeofday(&end,0);

    c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
    std::cout<<"Computation time set Solver: "<<c_time<<" s"<<std::endl;

    // nlp_solver.evaluate();

    // std::cout << std::setw(30) << "Objective: " << nlp_solver.output("f").getDescription() << std::endl;
    // std::cout << std::setw(30) << "Primal solution: " << nlp_solver.output("x").getDescription() << std::endl; 
  }
  //-------------------------------------------------------------------
  void SearchZones::extractPhList(std::vector<Eigen::MatrixXd*>& Ph_list)const
  {
    //Ph_list stores, for each of the h facets of the prototype GWS's convex hull, the vertices spanning the h-th facet

    uint H=grasp_->getGWS()->getNumFacets();
    uint K=grasp_->getGWS()->getDimension();

    Ph_list.resize(H);

    facetT* curr_f=grasp_->getGWS()->getConvexHull()->beginFacet().getFacetT();
    for (uint h=0; h<H; h++)
      {
	setT* vertices=curr_f->vertices;
	uint V=(uint)qh_setsize(vertices);
	Eigen::MatrixXd* Ph(new  Eigen::MatrixXd(K,V));

	//put the vertices of the current facet in the columns of the matrix Ph
	for(uint v=0;v<V ;v++)
	  {
	    coordT* p =((vertexT*)vertices->e[v].p)->point;
            for (uint k=0;k<K;k++) //Really shouldn't be doint that in a loop...
	      (*Ph)(k,v)=p[k];
	  }

	Ph_list[h]=Ph;
	curr_f=curr_f->next;
      }
  }
  //-------------------------------------------------------------------
  void SearchZones::mapFacetToFinger(const orgQhull::QhullFacet& facet,IndexList & finger_ids)const
  {
    finger_ids.clear();

    for (uint v=0; v<facet.vertices().size(); v++)
      {
	uint v_id=(uint)facet.vertices()[v].point().id();

	//determine to which finger the current vertex belongs
	for(uint n=0; n < grasp_->getNumFingers(); n++)
	  if(v_id <= map_vertex2finger_(n))
	    {
	      finger_ids.push_back(n);
	      break;
	    }
      }
    finger_ids.sort();
    finger_ids.unique();
  }
  //-------------------------------------------------------------------
  void SearchZones::extractWhiList(std::vector< ContactRegion * > const & conditioning_patches,std::vector<std::vector<Eigen::MatrixXd*> >& Wh_i_list)const
  {
    //Wh_i_list stores, for each of the h facets of the prototype GWS's convex hull, the wrench cones of the conditioning patches belonging to those fingers of the prototype grasp whose wrenches span the h-th facet

    uint H=grasp_->getGWS()->num_facets_;
    Wh_i_list.resize(H);
    IndexList finger_ids;
    orgQhull::QhullFacet facet=grasp_->getGWS()->conv_hull_.beginFacet();

    for (uint h=0; h<H; h++)
      {
	mapFacetToFinger(facet,finger_ids);
	std::vector<Eigen::MatrixXd*> Wh_i;
	for (IndexList::const_iterator f_it = finger_ids.begin(), end = finger_ids.end(); f_it != end; ++f_it)
	  for (uint j=0; j<conditioning_patches[*f_it]->size();j++)
	    {
	      Patch* patch=conditioning_patches[*f_it]->at(j);
	      for (IndexList::const_iterator p_it = patch->patch_ids_.begin(), end = patch->patch_ids_.end(); p_it != end; ++p_it)
		Wh_i.push_back(new Eigen::MatrixXd(*(grasp_->getFinger(*f_it)->getOWS()->getWrenchCone(*p_it)->getWrenches())));

	    }

	Wh_i_list[h]=Wh_i;
	facet=facet.next();
      }
  }
  //-------------------------------------------------------------------
  void SearchZones::computeConditionedHyperplanes(std::vector< ContactRegion * > const & conditioning_patches)
  {
    assert(grasp_->getGWS()->convHullComputed());
    uint H=grasp_->getGWS()->num_facets_;
    std::vector<std::vector<Eigen::MatrixXd* > > Wh_i_list(H);
    std::vector<Eigen::MatrixXd*> Ph_list(H);
    CasADi::IpoptSolver nlp_solver;


    struct timeval start, end;
    double c_time;

    if(tws_->getWrenchSpaceType() == Spherical)
      {
	std::cout<<"Error in SearchZones::computeConditionedHyperplanes() - Spherical Task Wrench Space based search zone computation not implemented yet!"<<std::endl;
      }
    else if(tws_->getWrenchSpaceType() == Discrete)
      {
	std::vector<Eigen::MatrixXd*> Ph_list;
        extractPhList(Ph_list); //ALLOCATES MEMORY!!! Should be done with shared pointers ...
        extractWhiList(conditioning_patches,Wh_i_list);//ALLOCATES MEMORY!!! Should be done with shared pointers ...

	createDiscreteTWSNLPSolver(Ph_list,Wh_i_list,nlp_solver);
	std::cout<<"TRYING TO SOLVE..."<<std::endl;
	gettimeofday(&start,0);
	nlp_solver.evaluate();

	// std::cout << std::setw(30) << "Objective: " << nlp_solver.output("f").getDescription() << std::endl;
	// std::cout << std::setw(30) << "Primal solution: " << nlp_solver.output("x").getDescription() << std::endl;

	gettimeofday(&end,0);
	c_time = end.tv_sec - start.tv_sec + 0.000001 * (end.tv_usec - start.tv_usec);
	std::cout<<"Computation time SOLVER: "<<c_time<<" s"<<std::endl;
	setTWSNLPHyperplanes(nlp_solver,Wh_i_list);
      }
    else
      {
	std::cout<<"Error in SearchZones::computeConditionedHyperplanes() - invalid Task Wrench Space type!"<<std::endl;
	return;
      }

    //clean up
    for (uint h=0; h<H; h++)
      {
	for (uint i=0; i<Wh_i_list[h].size();i++)
	  delete Wh_i_list[h][i];  

        delete Ph_list[h];
      }
  }
  //-------------------------------------------------------------------
  void SearchZones::addShiftedPrimitiveSearchZone(uint finger_id,vertexT const* curr_vtx)
  {
    setT *neighbor_facets=curr_vtx->neighbors;
    PrimitiveSearchZone* p_search_zone = new PrimitiveSearchZone();
    p_search_zone->hyperplane_ids_.reserve(qh_setsize(neighbor_facets));

    for(uint i=0; i < (uint)qh_setsize(neighbor_facets); i++)
      p_search_zone->hyperplane_ids_.push_back(((facetT*)neighbor_facets->e[i].p)->id);

    search_zones_[finger_id]->push_back(p_search_zone);
  }
  //-------------------------------------------------------------------
  void SearchZones::clear()
  {
    for(uint finger_id=0;finger_id < num_search_zones_;finger_id++)
      {
	for(uint prim_sz_id=0; prim_sz_id < search_zones_[finger_id]->size(); prim_sz_id++)
	  delete (*search_zones_[finger_id])[prim_sz_id];            

	delete search_zones_[finger_id];
      }
    search_zones_.clear();
    num_search_zones_ = 0;
    search_zones_computed_ = false;

  }
  //-------------------------------------------------------------------
  void SearchZones::initializeSearchZones()
  {
    if(search_zones_computed_) //Clear up possible previously computed search zones
      clear();

    num_search_zones_=grasp_->getNumFingers();
    search_zones_.reserve(num_search_zones_);
    map_vertex2finger_.resize(1,num_search_zones_);

    for(uint finger_id=0; finger_id < num_search_zones_; finger_id++)
      {
	search_zones_.push_back(new SearchZone());
	//reserve for the maximum number of primitive search zones for each finger (i.e all wrenches are vertexes of the gws and span a primitive search zone)
	search_zones_[finger_id]->reserve(grasp_->getFinger(finger_id)->getCenterPointPatch()->patch_ids_.size()
					  *grasp_->getFinger(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches());

	//the elements of map_vertex2finger_ describe the relation between a vertex-index from the grasp wrench space to the finger this vertex belongs.
	//I.e. if the value of vertex id <= than map_vertex2finger_(0), then the vertex is a wrench belonging to the first finger, if vertex id <= than 
	//map_vertex2finger_(1) the vertex belongs to the second finger and so on. E.g, a 3-fingered grasp comprising patches with 2 contact points each and 10 
	//primitive wrenches in each contact point would have map_vertex2finger_=[19 39 59] 
	map_vertex2finger_(finger_id)=grasp_->getFinger(finger_id)->getCenterPointPatch()->patch_ids_.size()
	  *grasp_->getFinger(finger_id)->getContactModel()->getLimitSurface()->getNumPrimitiveWrenches()-1;
	if (finger_id > 0)
	  map_vertex2finger_(finger_id)=map_vertex2finger_(finger_id)+map_vertex2finger_(finger_id-1)+1;
      }
  }
  //--------------------------------------------------------------------
  void SearchZones::computeShiftedSearchZones()
  {
    //Make sure the TWS is valid
    assert(tws_.get() != NULL);
    assert(tws_->getDimension() == grasp_->getGWS()->getDimension());
    assert(grasp_->getGWS()->getOcInsphereRadius() >= tws_->getOcInsphereRadius());

    //Make sure the grasp is valid
    assert(grasp_->getGWS()->containsOrigin());

    initializeSearchZones();
    computeShiftedHyperplanes();
    computePrimitiveSearchZones();
 
    search_zones_computed_=true;
  }
  //--------------------------------------------------------------------
  void SearchZones::computeConditionedSearchZones(std::vector< ContactRegion * > const & conditioning_patches)
  {
    //Make sure the TWS is valid
    assert(tws_.get() != NULL);
    assert(tws_->getDimension() == grasp_->getGWS()->getDimension());
    assert(grasp_->getGWS()->getOcInsphereRadius() >= tws_->getOcInsphereRadius());

    //Make sure the grasp is valid
    assert(grasp_->getGWS()->containsOrigin());

    //Make sure the conditioning patches are valid
    assert(conditioning_patches.size()==grasp_->getNumFingers());

    initializeSearchZones();
    computeConditionedHyperplanes(conditioning_patches);
    computePrimitiveSearchZones();
 
    search_zones_computed_=true;
  }
  //--------------------------------------------------------------------
  void SearchZones::computePrimitiveSearchZones()
  { 
    //iterate through the vertices of the GWS and create the appropriate Primitive Search Zone for each 
    //vertex and push it back on the search zone of its corresponding finger

    orgQhull::QhullVertex curr_vtx=grasp_->getGWS()->conv_hull_.beginVertex();

    for (uint vtx_count=0; vtx_count < grasp_->getGWS()->num_vtx_;vtx_count++)
      {
	//determine to which finger the current vertex belongs and create an appropriate primitive search zone and push it on the list for the corresponding
	//finger
	for(uint finger_id=0; finger_id < num_search_zones_; finger_id++)
	  {
	    if((uint)curr_vtx.point().id() <= map_vertex2finger_(finger_id))
	      {
		addShiftedPrimitiveSearchZone(finger_id,curr_vtx.getVertexT());
		break;
	      }
	  }

	curr_vtx=curr_vtx.next();
      }
  }
  //--------------------------------------------------------------------
  Eigen::Matrix<double,Eigen::Dynamic,6> const* SearchZones::getHyperplaneNormals()const{return &hyperplane_normals_;}
  //--------------------------------------------------------------------
  Eigen::VectorXd const* SearchZones::getHyperplaneOffsets()const{return &hyperplane_offsets_;}
  //--------------------------------------------------------------------
  const GraspPtr SearchZones::getGrasp()const{return grasp_;}
  //--------------------------------------------------------------------
  WrenchSpacePtr SearchZones::getTaskWrenchSpace()const
  {
    assert(search_zones_computed_);
    return tws_;
  }
  //--------------------------------------------------------------------
  void SearchZones::resetPrimitiveSearchZones(uint sz_id)
  {
    for(uint psz_id=0; psz_id < search_zones_[sz_id]->size(); psz_id++)
      (*search_zones_[sz_id])[psz_id]->satisfied_wc_ids_.resize(0);  
  }
  //--------------------------------------------------------------------
  void SearchZones::setTaskWrenchSpace(WrenchSpacePtr tws)
  {
    //make sure the TWS is valid
    assert(tws.get() != NULL);

    tws_=tws;
  }
  //--------------------------------------------------------------------
  void SearchZones::resetSearchZones()
  {
    for(uint sz_id=0;sz_id < num_search_zones_; sz_id++)
      resetPrimitiveSearchZones(sz_id);
  }
  //--------------------------------------------------------------------
  //--------------------------------------------------------------------
}//namespace ICR
