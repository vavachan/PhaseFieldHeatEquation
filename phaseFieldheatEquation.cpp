/* These are general */
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <random>
#include <cuda_runtime.h>   // for cudaMemGetInfo profiling probes


/* openLB header files */ 
#include "olb3D.h"
#include "olb3D.hh" 

#include<time.h>

//namespace olb::descriptors{
//	struct PREV_PHI		    : public FIELD_BASE<1,  0, 0> { };
//}
using namespace olb; 
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::graphics;
using namespace olb::dynamics;

//typedef Kernel::Point_3 Point_3;
//typedef Kernel::FT FT;
using T = FLOATING_POINT_TYPE;

using DESCRIPTOR =  D3Q15<POROSITY,PREV_PHI>;
using TDESCRIPTOR = D3Q7<>;


//using BulkDynamics = olb::dynamics::Tuple<
//	T, DESCRIPTOR,
//	momenta::BulkTuple,
//	equilibria::SecondOrder,
//	collision::BGK,
//	forcing::Guo<momenta::Forced>
//	>;



const T  rho    = 1000;     //Kg/m^3 // We will keep rho same for water and ice

const T  C_p_W    = 4184.0;   // m^2 K^-1 s^-2
const T  k_cond_W = 0.58;     // Kg⋅m⋅s^−3⋅K^−1

const T  C_p_I    = C_p_W;//2108.0;
const T  k_cond_I = 2.2; 

const T  L_heat = 334000.0; // m^2.s^-2 
const T  T_m    = 273.15;   // K

const T  k_cond_S = 2.2; 

const T  kappa  = k_cond_I / (rho * C_p_W); // Kg⋅m⋅s^−3⋅K^−1 / ( Kg.m^-3.m^2 K^-1 s^-2) = m^2/s

const T interfaceFreeEnergyIceWater = 3*1e-2 ; // J / m^2
const T heatCapacityVolumetric = C_p_W*910; // J / [m^3 . K]
const T latentHeatVolumetric = L_heat*910;  // J / m^3
const T thermalCapilaryLength = 4.*interfaceFreeEnergyIceWater*heatCapacityVolumetric/(latentHeatVolumetric*latentHeatVolumetric)*T_m; // [gamma.C_p / L^2] * T_m : [m]


T Temp;// = 272.85;

T Delta;// = abs(C_p_W * (Temp - T_m) / L_heat);

const T a1 = 0.88388347648;
const T a2 = 0.6267;

const T lengthX = 0.15*1e-3;
const T lengthY = 0.15*1e-3;
const T lengthZ = 0.15*1e-3;

int N = 450;

const T dx= lengthX/N;
//const T dt= 1.;

const T W = 3*dx;
const T W_sq = W * W;         // W^2 (for Eq 1)
const T lambda =  6;
const T tau = lambda*(a2*W_sq)/kappa; 

const T eta_theta= 1.;//descriptors::invCs2<T,TDESCRIPTOR>()*kappa*(dt/(dx*dx))+1./2.;
const T dt = (eta_theta-1./2.)*(dx*dx)/(descriptors::invCs2<T,TDESCRIPTOR>()*kappa);
const T eta_phi = descriptors::invCs2<T,DESCRIPTOR>()*W_sq/tau*(dt/(dx*dx))+1./2.;




/// Logs current GPU memory usage and the delta since the last call.
/// Call this right after each major construction step to see which
/// stage of setup/run is actually consuming memory.
void reportGpuMem(const std::string& label)
{
	static size_t prevUsedB = 0;
	size_t freeB = 0, totalB = 0;
	cudaMemGetInfo(&freeB, &totalB);
	size_t usedB = totalB - freeB;
	if (singleton::mpi().getRank() == 0) {
		const double GiB = 1024.0*1024.0*1024.0;
		std::cout << "[GPU MEM] " << label << ": "
		          << (usedB/GiB) << " GiB used"
		          << "  (delta " << ((double)usedB - (double)prevUsedB)/GiB << " GiB)"
		          << "  / " << (totalB/GiB) << " GiB total"
		          << std::endl;
	}
	prevUsedB = usedB;
}

template <typename T, typename S, typename DESCRIPTOR>
class TanhSphere : public AnalyticalF3D<T, S>
{
	private:
		T *x,*y,*z,radius, ls; // ls is the width of the border . 
		int nNucli;
		//UnitConverter<T, DESCRIPTOR> const &converter;
		//int iT, tempCase;

	public:
		TanhSphere(T *x_, T *y_,T *z_, T radius_, T ls_, int nNucli_) : AnalyticalF3D<T, S>(3),
		x(x_), y(y_), z(z_), radius(radius_), ls(ls_), nNucli(nNucli_)
	{
		this->getName() = "TanhSphere";
	};

		bool operator()(T output[1], const S X[3])
		{
			//T res = converter.getResolution();
			//T *outputs = new T[nNucli];
			std::vector<T> outputs(nNucli);         // ← stack-managed, no leak

			for(int i=0; i<nNucli; i++)
			{
				T real_distx = X[0];
				T real_disty = X[1];
				T real_distz = X[2];
				T distance = sqrt((real_distz-z[i])*(real_distz-z[i])+(real_distx-x[i])*(real_distx-x[i])+(real_disty-y[i])*(real_disty-y[i]));
				T argument = (radius - distance)/(sqrt(ls));

				outputs[i] = util::tanh(argument);

			}
			//output[0] = *std::max_element(outputs,outputs+nNucli);
			output[0] = *std::max_element(outputs.begin(), outputs.end());


			return true;
		};
};

template <unsigned D>
const auto periodicity_func = [](){
	if constexpr (D==3) {
		return Vector<bool,3>(true,true,false);
	} else {
		return Vector<bool,2>(false,false);
	}
};


Vector<T,2> computeAveragePorosity(SuperLattice<T,DESCRIPTOR>& sLattice)
{
    int  localSumSolid = 0;
    int  localSumLiquid = 0;
    int  localSum = 0;
    int  localCount  = 0;

    for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
        auto& block = sLattice.getBlock(iC);
        block.forCoreSpatialLocations([&](LatticeR<3> loc) {
            auto cell = block.get(loc);
            T por = cell.template getField<descriptors::POROSITY>();

            if (por >= T{0.99}) 
	    {
                T rho{};
                T u[3]{};                    // raw array, not Vector
                cell.computeRhoU(rho, u);
		if(rho > 0.0)
		{
			localSumSolid += 1;
		}
		else
		{
			localSumLiquid += 1;
		}
                localCount++;
            }
        });
    }
    // reduceAndBcast reduces in-place on each rank
    singleton::mpi().reduceAndBcast(localSumLiquid, MPI_SUM);
    singleton::mpi().reduceAndBcast(localSumSolid, MPI_SUM);
    singleton::mpi().reduceAndBcast(localCount,  MPI_SUM);

    return Vector<T,2>{localSumSolid*1./localCount,
                       localSumLiquid*1./localCount};
}

void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,SuperGeometry<T,3>& superGeometry, bool restart)
{
    T L = converter.getPhysDeltaX();

    OstreamManager clout(std::cout,"Prepare Geometry");
    superGeometry.rename(0, 1);

    Vector<T,3> extend_Bottom(lengthX, lengthY, L);
    Vector<T,3> origin_Bottom(0., 0., -L); 
    IndicatorCuboid3D<T> bottombox(extend_Bottom, origin_Bottom);

    Vector<T,3> extend_Top(lengthX, lengthY, L);
    Vector<T,3> origin_Top(0., 0., lengthZ );
    IndicatorCuboid3D<T> topbox(extend_Top, origin_Top);
    
    superGeometry.rename(1, 2, bottombox);

    superGeometry.rename(1, 3, topbox);

    superGeometry.clean(); 

    superGeometry.checkForErrors();
    superGeometry.print();
    
    clout << "Geometry setup complete: Bulk=1, Bottom=2, Top=3" << std::endl;
}
void prepareLattices(SuperGeometry<T,3>& superGeometry,
		SuperLattice<T,DESCRIPTOR>& phaseField,
		SuperLattice<T,TDESCRIPTOR>& heatEquation
//		SuperLattice<T,DESCRIPTOR>& fluid,
//		SuperLattice<T,TDESCRIPTOR>& SS,
		//UnitConverter<T,DESCRIPTOR> const& converter
		)
{


	OstreamManager clout(std::cout,"prepare Lattice");

	// Both lattices evolve through the diffusion dynamics.
	// ALL DYNAMICS IS DEFINED HERE
	//const T omega = converter.getLatticeRelaxationFrequency();


	/* PHASE FIELD */
	phaseField.defineDynamics<DiffusionDynamics>(superGeometry,1);
	//boundary::set<boundary::LocalVelocity>(phaseField, superGeometry, 2);
	//boundary::set<boundary::LocalPressure>(phaseField, superGeometry, 3);
	phaseField.defineDynamics<BounceBack>(superGeometry, 2);
	phaseField.defineDynamics<BounceBack>(superGeometry, 3);
	/*########################*/

	/* HEAT EQUATION */
	heatEquation.defineDynamics<DiffusionDynamics>(superGeometry,1);
	//boundary::set<boundary::AdvectionDiffusionDirichlet>(heatEquation, superGeometry, 2);
	//boundary::set<boundary::AdvectionDiffusionDirichlet>(heatEquation, superGeometry, 3);
	boundary::set<boundary::LocalPressure>(heatEquation, superGeometry, 2);
	boundary::set<boundary::LocalPressure>(heatEquation, superGeometry, 3);

	//setLocalPressureBoundary(heatEquation, omega, superGeometry, 2);
	//setLocalPressureBoundary(heatEquation, omega, superGeometry, 3);
	//heatEquation.defineDynamics<BounceBack>(superGeometry,2);
	//heatEquation.defineDynamics<BounceBack>(superGeometry,3);
	/*########################*/

	/*FLUID*/
	//fluid.defineDynamics<BulkDynamics>(superGeometry,1);
	//boundary::set<boundary::LocalVelocity>(fluid, superGeometry, 2);
	//boundary::set<boundary::LocalPressure>(fluid, superGeometry, 3);

	//setLocalVelocityBoundary(fluid, omega, superGeometry, 2);
	//setLocalPressureBoundary(fluid, omega, superGeometry, 3);
	/*********************/

	/*Super Saturation Equation*/
	//SS.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry,1);
	//boundary::set<boundary::AdvectionDiffusionDirichlet>(SS, superGeometry.getMaterialIndicator({2, 3}));
	/**************************/



	clout<<"Lattice preparation done"<<std::endl;

}

void initializeLattice(SuperGeometry<T,3>& superGeometry,
		SuperLattice<T,DESCRIPTOR>& phaseField,
		SuperLattice<T,TDESCRIPTOR>& heatEquation,
		//SuperLattice<T,DESCRIPTOR>& fluid,
		//SuperLattice<T,TDESCRIPTOR>& SS,
		T *x, T *y, T *z,
		int nNucli)
		//UnitConverter<T,DESCRIPTOR> const& converterFLUID)
{
	OstreamManager clout(std::cout,"initialize Lattice");
	clout<<"Initialization of Lattices"<<std::endl;
	auto bulkIndicator = superGeometry.getMaterialIndicator({1,2,3});
	auto bulkIndicatorHE = superGeometry.getMaterialIndicator({1});
	auto top = superGeometry.getMaterialIndicator({3});
	auto bottom = superGeometry.getMaterialIndicator({2});

	TanhSphere <T,T,DESCRIPTOR> tanhs(x,y,z,30*dx,2*W*W,nNucli);

	std::vector<T> velocity(3,T(0.));   // Velocity is zero
	AnalyticalConst3D<T,T> uF(velocity);
	//AnalyticalConst3D<T,T> temp_HI(1-Delta);
	//AnalyticalConst3D<T,T> temp_LO(1-Delta);
	AnalyticalConst3D<T,T> underCooling(1-Delta);
	AnalyticalConst3D<T,T> phiLiquid(-1.);
	AnalyticalConst3D<T,T> phiSolid(1.);

	/* HEAT EQUATION */

	heatEquation.defineRhoU(bulkIndicator,underCooling,uF);
	heatEquation.iniEquilibrium(bulkIndicator,underCooling,uF);

	heatEquation.defineRhoU(bottom,underCooling,uF);
	heatEquation.iniEquilibrium(bottom,underCooling,uF);

	heatEquation.defineRhoU(top,underCooling,uF);
	heatEquation.iniEquilibrium(top,underCooling,uF);

	clout<<"Initialized heat equation"<<std::endl;
	/*#######################*/

	phaseField.defineRhoU(bulkIndicator,phiLiquid,uF);
	phaseField.iniEquilibrium(bulkIndicator,phiLiquid,uF);
	//phaseField.defineRhoU(bulkIndicator,phiSolid,uF);
	//phaseField.iniEquilibrium(bulkIndicator,phiSolid,uF);

	//phaseField.defineRhoU(bulkIndicator,tanhs,uF);
	//phaseField.iniEquilibrium(bulkIndicator,tanhs,uF);

	phaseField.defineRhoU(top, phiSolid, uF);
	phaseField.iniEquilibrium(top, phiSolid, uF);
	//phaseField.defineRhoU(bottom, phiSolid, uF);
	//phaseField.iniEquilibrium(bottom, phiSolid, uF);

	phaseField.initialize();
	heatEquation.initialize();

	clout<<"Lattice initialization done"<<std::endl;
}



int readConfigToXYZ(T *&X, T *&Y, T *&Z, T *&RAD, std::string filename)
{
	OstreamManager clout(std::cout,"readConfigToXYZ");
	clout<<"starting ..."<<std::endl;
	std::ifstream config(filename);
	std::string STR;
	int nAtoms = 15000;
	T scale=1e-3;
//	for(int i=0; i <9;i++)
//	{
//		std::getline(config,STR);
//		if(i==3)
//			nAtoms = std::stoi(STR);
//		//std::cout<<STR<<"\n";
//	}

	X = new T[nAtoms];
	Y = new T[nAtoms];
	Z = new T[nAtoms];
	RAD = new T[nAtoms];
	//radius = 1*sigma/deltax;
	for(int i=0; i < nAtoms; i++)
	{
		config>>X[i]>>Y[i]>>Z[i]>>RAD[i];
		X[i] = X[i]*scale;
		Y[i] = Y[i]*scale;
		Z[i] = Z[i]*scale;
		RAD[i] = RAD[i]*scale;
		//clout<<i<<"\t"<<X[i]<<"\t"<<Y[i]<<"\t"<<Z[i]<<RAD[i]<<std::endl;
	}
	config.close();
	clout<<"done..."<<std::endl;
	return nAtoms;
}
void findNNucliLocations(T *X, T *Y, T *Z, T *RAD, int nAtoms,
	       	T *x, T *y, T *z, int nNucli)
{
	OstreamManager clout(std::cout,"findNNucliLocations");
	clout<<"starting ..."<<std::endl;
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	const int N = nNucli;         
	const T Lx = lengthX;  
	const T Ly = lengthY;   
	const T Lz = lengthZ;   

	unsigned int common_seed = 12345;

	std::mt19937 gen(common_seed);

	std::uniform_real_distribution<T> dist_x(0.0, Lx);
	std::uniform_real_distribution<T> dist_y(0.0, Ly);
	std::uniform_real_distribution<T> dist_z(0.0, Lz);

	for (int i = 0; i < N; ++i) {
		x[i] = dist_x(gen);
		y[i] = dist_y(gen);
		z[i] = dist_z(gen);
		clout<<x[i]<<"\t"<<y[i]<<"\t"<<z[i]<<std::endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	clout<<"done ..."<<std::endl;
}

int main( int argc, char* argv[] )
{
	//olbInit( &argc, &argv );
	olb::initialize( &argc, &argv );

	srand(static_cast<unsigned int>(time(0)));
	singleton::directories().setOutputDir( argv[1]);
	OstreamManager clout(std::cout,"main");
	Temp = std::stof(argv[2]);
	//int rand (void);
	clout<<"Simulation start.."<<std::endl;

	clout<<"OpenLB parameters"<<std::endl;

	UnitConverterFromResolutionAndRelaxationTime<T, TDESCRIPTOR> const converterHE(
			int {1},                        // resolution: number of voxels per charPhysL
			(T)   eta_theta,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
			(T)   dx*1,       // charPhysLength: reference length of simulation geometry
			(T)   0.02*(dx/dt),                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
			(T)   kappa, // physViscosity: physical kinematic viscosity in __m^2 / s__
			(T)   1.0                       // physDensity: physical density in __kg / m^3__
			);
	UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterPF(
			int {1},                        // resolution: number of voxels per charPhysL
			(T)   eta_phi,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
			(T)   dx*1,       // charPhysLength: reference length of simulation geometry
			(T)   0.02*(dx/dt),                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
			(T)   W_sq/tau, // physViscosity: physical kinematic viscosity in __m^2 / s__
			(T)   1.0                       // physDensity: physical density in __kg / m^3__
			);
	//const UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converterFLUID(
	//		int{1},                      // resolution: number of voxels per charPhysL
	//		(T)   eta_fluid,
	//		(T)   dx*1,                // charPhysLength: reference length of simulation geometry
	//		(T)   0.4,                // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
	//		(T)   viscositySim,             // physViscosity: physical kinematic viscosity in __m^2 / s__
	//		(T)   1000		// physDensity: physical density in __kg / m^3__
	//		);
	//UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterSS(
	//		int {1},                        // resolution: number of voxels per charPhysL
	//		(T)   eta_U,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
	//		(T)   dx*1,       // charPhysLength: reference length of simulation geometry
	//		(T)   0.02*(dx/dt),                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
	//		(T)   D_sim, // physViscosity: physical kinematic viscosity in __m^2 / s__
	//		(T)   1.0                       // physDensity: physical density in __kg / m^3__
	//		);

	clout<<eta_phi<<"\t"<<eta_theta<<std::endl;
	//return 0;
	converterHE.print();
	converterPF.print();
	//const T omega_FLUID = converterFLUID.getLatticeRelaxationFrequency();
	//const T omega_SS= converterSS.getLatticeRelaxationFrequency();
	clout << "Thermal Capilary length (actual) [d_0] \t"<< thermalCapilaryLength << " m"<<std::endl;
	clout << "Thermal Capilary length (numerical) [d_num] \t"<< a1*W/lambda<< " m"<<std::endl;

	clout << "Temp we WANT to simulate in \t"<< Temp <<"K"<<std::endl;
	
	T Delta_Real  = abs(C_p_W * (Temp - T_m) / L_heat);
	T Delta_numerical = ((a1*W/lambda)/thermalCapilaryLength)*Delta_Real;
	T d_num_d_0_ratio = (a1*W/lambda)/thermalCapilaryLength;
	T dTemp_numerical = L_heat/C_p_W * Delta_numerical;

	clout << "Temp we SHOULD simulate in \t"<<  T_m - dTemp_numerical <<"K"<< std::endl;

	Delta = Delta_numerical;
	clout << "critical nucleation radius (numerical) \t"<<1.768*W/(lambda*Delta)/dx<<" l.u; \t"<<1.768*W/(lambda*Delta)<<" m" <<std::endl;


	T L=converterHE.getPhysDeltaX();
	bool restart = false;

	// These two define the whole region
	Vector<T,3> extend(lengthX,lengthY,lengthZ);
	Vector<T,3> origin;//(-5,-2.5,-2.5);
	IndicatorCuboid3D<T> region(extend,origin);

#ifdef PARALLEL_MODE_MPI
	const int noOfCuboids = singleton::mpi().getSize();
#else
	const int noOfCuboids = 1;
#endif

	CuboidDecomposition<T,3> regionGeometry(region,L,noOfCuboids); // Define the whole region
	regionGeometry.setPeriodicity({true,true,false});  // Periodic boundary condition in x and y directions 

	HeuristicLoadBalancer<T> loadBalancer( regionGeometry );
	SuperGeometry<T,3> superGeometry(regionGeometry,loadBalancer);
	prepareGeometry(converterPF,superGeometry,restart);

	clout<<"geometry created"<<std::endl;
	reportGpuMem("after geometry");
	/* One lattice will solve the phase field equation 
	 * and the another one will solve the heat equation.
	 * 
	 * The equation are on the same geometry, and they also 
	 * follow the same dynamics essentially. 
	 *
	 * The coupling will add the corresponding terms through 
	 * which the lattices interact with each other. 
	 */

	SuperLattice<T,DESCRIPTOR> phaseField(superGeometry);
	SuperLattice<T,TDESCRIPTOR> heatEquation(superGeometry);



	prepareLattices(superGeometry,phaseField,heatEquation);//,fluid,SS,converterFLUID);

	phaseField.setParameter<descriptors::OMEGA>(1./eta_phi);//omega_PF);
	heatEquation.setParameter<descriptors::OMEGA>(1./eta_theta);
	T *X=nullptr;
	T *Y=nullptr;
	T *Z=nullptr;
	T *RAD=nullptr;
	int nAtoms = readConfigToXYZ(X,Y,Z,RAD,"N15000LogNormalDistributionJUSTPOSITIONS.dat");


	int nNucli = 1;
	T *x = new T[nNucli];
	T *y = new T[nNucli];
	T *z = new T[nNucli];
	//findNNucliLocations(X,Y,Z,RAD,nAtoms,x,y,z,nNucli);
	x[0] = lengthX/2.;
	y[0] = lengthY/2.;
	z[0] = lengthZ/2.;
	initializeLattice(
			superGeometry,
			phaseField,
			heatEquation,
			x,y,z,
			nNucli);


	reportGpuMem("after initialization");


	SuperPorousElementEmbeddingO fsiEmbeddingO(phaseField);
	reportGpuMem("before registration");
	fsiEmbeddingO.registerElementType<SpherePorosityF>(superGeometry.getMaterialIndicator(1));
	reportGpuMem("before sphere");
  	for (unsigned iSphere=0; iSphere < nAtoms; ++iSphere) 
  	{
  		Vector<T,3> spherePhysR ={X[iSphere],Y[iSphere],Z[iSphere]};
  		if(iSphere%(1000)==0)
  		{
  			clout<<"Added "<<iSphere<<" spheres"<<std::endl;
			reportGpuMem("after +thousand");
  		}
  		auto sphereParams = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<3>>(
  				fields::fsi::ELEMENT_TAG{},   1+iSphere,
  				fields::fsi::ELEMENT_PIVOT{}, spherePhysR,
  				SpherePorosityF::RADIUS{},    RAD[iSphere]
  				);
  		fsiEmbeddingO.add(sphereParams);
  	}



	phaseField.setParameter<fields::converter::PHYS_DELTA_X>(converterPF.getPhysDeltaX());

	phaseField.setProcessingContext(ProcessingContext::Simulation); // ← and this
	heatEquation.setProcessingContext(ProcessingContext::Simulation); // ← and this




	reportGpuMem("before embedding");

	fsiEmbeddingO.initialize();
	fsiEmbeddingO.apply();
	phaseField.communicate();

	reportGpuMem("after embedding");
	{
		auto& c = phaseField.getCommunicator(stage::Full{});
		c.template requestField<descriptors::POROSITY>();
		c.template requestField<fields::fsi::ELEMENT_TAG>();
		c.exchangeRequests();
	}


	clout<<"done adding sphere"<<std::endl;
	fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_PIVOT>>(ProcessingContext::Simulation);

	//SuperLatticeVelocity3D <T,DESCRIPTOR> velocityPF(phaseField);
	SuperLatticeDensity3D <T,DESCRIPTOR> densityPF(phaseField);
	//SuperLatticeField3D <T,DESCRIPTOR,POPULATION> populationPF(phaseField);
	SuperLatticeField3D <T,DESCRIPTOR,POROSITY> externalPor(phaseField);//,"porosity");


	//SuperLatticeVelocity3D <T,DESCRIPTOR> velocityHE(heatEquation);
	SuperLatticeDensity3D <T,TDESCRIPTOR> densityHE(heatEquation);
	//SuperLatticeField3D <T,DESCRIPTOR,POPULATION> populationHE(heatEquation);

	SuperVTMwriter3D<T> vtmWriterPF("PhaseField");
	//SuperVTMwriter3D<T> vtmWriterHE("HeatEquation");


	vtmWriterPF.createMasterFile();
	//vtmWriterHE.createMasterFile();


	vtmWriterPF.addFunctor(densityPF);
	//vtmWriterPF.addFunctor(velocityPF);
	vtmWriterPF.addFunctor(externalPor);
	vtmWriterPF.addFunctor(densityHE,"heatEquation");


	//vtmWriterHE.addFunctor(densityHE);
	//vtmWriterHE.addFunctor(velocityHE);
	//vtmWriterHE.addFunctor(populationHE);


	vtmWriterPF.write(-2);
	SuperLatticeCoupling couplingFirst(PhaseFieldHeatEquationCouplingFirstPart<T>{},
			names::A{},phaseField,
			names::B{},heatEquation
			);

	SuperLatticeCoupling couplingSecond(PhaseFieldHeatEquationCouplingSecondPart<T>{},
			names::A{},phaseField,
			names::B{},heatEquation
			);


	couplingFirst.setParameter<PhaseFieldHeatEquationCouplingFirstPart<T>::TAU>(tau);
	couplingFirst.setParameter<PhaseFieldHeatEquationCouplingFirstPart<T>::DELTA_T>(dt);
	couplingFirst.setParameter<PhaseFieldHeatEquationCouplingFirstPart<T>::LAMBDA>(lambda);
	couplingFirst.setParameter<PhaseFieldHeatEquationCouplingFirstPart<T>::OMEGAPFE>(1./eta_phi);

	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::thermalConductivityWater>(k_cond_W);
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::thermalConductivityIce>(k_cond_I);
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::thermalConductivitySoil>(k_cond_S);
                                                                           
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::heatCapacityWater>(C_p_W);
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::heatCapacityIce>(C_p_I);
                                                                          
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::density>(rho);
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::etaPhi>(eta_phi);
                                                                         
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::deltax>(dx);
	couplingSecond.setParameter<PhaseFieldHeatEquationCouplingSecondPart<T>::deltat>(dt);


	//auto bulkIndicatorPF = superGeometry.getMaterialIndicator({2,3});
	auto topBottomIndicator= superGeometry.getMaterialIndicator({2,3});
	std::vector<T> velocity(3,T(0.));   // Velocity is zero
	AnalyticalConst3D<T,T> uF(velocity);


	T maxPhysT = 1000*60*1e-6;
	util::Timer<T> timer(converterPF.getLatticeTime(maxPhysT),
			superGeometry.getStatistics().getNvoxel());
	timer.start();
	int vtkWriteInter = converterHE.getLatticeTime(5*60*1e-6);
	int updateInter   = converterHE.getLatticeTime(5*60.*1e-6);
	int statInter     = converterHE.getLatticeTime(1*60.*1e-6);

	const std::string avgPorosityFileName = "averagePorosity.txt";
	std::ofstream avgPorosityOut(argv[1]+avgPorosityFileName, std::ios::out | std::ios::trunc);

	T min_temp = Delta ;
	T dT = min_temp; // below 273.15 K 
	T max_temp = T_m - 272.1;
	T DeltaTemp = (max_temp - min_temp)/(converterHE.getLatticeTime(maxPhysT/2.)/updateInter);

	clout << min_temp << "\t"<< max_temp <<"\t" << DeltaTemp << std::endl;
	AnalyticalConst3D<T,T> temp(1-Delta);
	heatEquation.defineRhoU(topBottomIndicator,temp,uF);
	heatEquation.iniEquilibrium(topBottomIndicator,temp,uF);
	heatEquation.template setProcessingContext<Array<momenta::FixedDensity::RHO>>(
			ProcessingContext::Simulation);
	clout<<Delta <<std::endl;
	for(int PF_time=0; PF_time*dt < maxPhysT; PF_time++)
	{
		if(PF_time % updateInter== 0 and PF_time > 0)
		{
			if(PF_time*dt > maxPhysT/2.)
				dT = dT - DeltaTemp;
			else
				dT = dT + DeltaTemp;
			clout<<"updating temp"<<std::endl;
			//clout<<PF_time << "\t"<< dT << std::endl;
			T delta = abs(C_p_I * (d_num_d_0_ratio * dT) / L_heat);
			AnalyticalConst3D<T,T> temp(1-Delta);
			heatEquation.defineRhoU(topBottomIndicator,temp,uF);
			heatEquation.iniEquilibrium(topBottomIndicator,temp,uF);
			heatEquation.template setProcessingContext<Array<momenta::FixedDensity::RHO>>(
					ProcessingContext::Simulation);

		}

		couplingFirst.execute();
		phaseField.collideAndStream();


		couplingSecond.execute();
		heatEquation.collideAndStream();

		if(PF_time == 0){
			reportGpuMem("after first collideAndStream step (steady-state footprint)");
		}

		if(PF_time%statInter==0)
		{
			timer.update(PF_time);
			timer.printStep();
			if(PF_time >0)
			{
				if (std::isnan(phaseField.getStatistics().getAverageRho())) {
					throw std::runtime_error("Simulation diverged.");
				}
			}
			phaseField.setProcessingContext(ProcessingContext::Evaluation);
        		//vtmWriterPF.write(PF_time);
			auto averagePorosity = computeAveragePorosity(phaseField);
			phaseField.setProcessingContext(ProcessingContext::Simulation); // ← and this
			if(singleton::mpi().getRank() == 0)
			{
				std::cout<< PF_time*dt <<"\t"<<averagePorosity<<"\t"<<dT<<"\n";
				avgPorosityOut << PF_time*dt <<"\t"<<averagePorosity[0]<<"\t"<<averagePorosity[1]<<"\t"<<dT<<std::endl;
			}
		}

		if(PF_time%vtkWriteInter==0)
		{
			phaseField.executePostProcessors(stage::Evaluation{});
			heatEquation.executePostProcessors(stage::Evaluation{});

			phaseField.setProcessingContext(ProcessingContext::Evaluation);
			heatEquation.setProcessingContext(ProcessingContext::Evaluation);

        		vtmWriterPF.write(PF_time);
        		//vtmWriterHE.write(PF_time);

			phaseField.save("phaseFieldEquation");
			heatEquation.save("heatEquation");


			phaseField.setProcessingContext(ProcessingContext::Simulation); 
			heatEquation.setProcessingContext(ProcessingContext::Simulation);

		}
	}

	// stop MPI environment
	MPI_Finalize();

	return 0;
} 
