#ifndef __MSU_BOLTZMANN_H__
#define __MSU_BOLTZMANN_H__

#include "msu_boltzmann/msupart.h"
#include "msu_commonutils/commondefs.h"
#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_commonutils/parametermap.h"
#include "msu_boltzmann/mutinfo.h"
#include "msu_commonutils/decay_nbody.h"
#include "msu_boltzmann/acceptance.h"
#include "msu_eos/resonances.h"
#include "msu_commonutils/log.h"
#include "msu_sampler/part.h"
#include "msu_boltzmann/balancearrays.h"
#include "msu_boltzmann/inelastic.h"
#include "msu_boltzmann/decay.h"

using namespace std;

class CpartList;
class CMSU_Decay;

class ChadronCount{
public:
	long long int Npi,NK,NB,NN,NLambda,NSigma,NXi,NOmega;
	ChadronCount(){
		Npi=NK=NB=NN=NLambda=NSigma=NXi=NOmega=0;
	}
};

class CMSU_Boltzmann{
public:
	int nevents;
	CparameterMap *parmap;
	CHBChargeMap chargemap;
	CMSUPartMap DeadPartMap;
	CMSUPartMap PartMap;		//!< A C++ map for active CMSUPart objects in the model.

	/*!
	This map is used to schedule and organize the various actions that the model must perform in time order. It contains all actions (as CAction objects) that have yet to occur, and the map's key is the boost-invariant time \f$\tau\f$ at which the action is scheduled to occur.
	*/
	CActionMap ActionMap;
	CActionMap DeadActionMap; // action objects not in line to be processed
	CresList *reslist;	//!< The CresList instance for the mucalodel (dynamically allocated).
	CInelasticList *inelasticlist;	//!< The CInelasicList instance for the model (dynamically allocated).
	
	int NXY;	//!< Determines size of mesh. The mesh size is \f$(2NXY,2NXY, 2NETA)\f$.
	int NETA;
	bool RESONANCE_DECAYS,BFCALC;
	int HYDRO_OCTANT_SYMMETRY;
	double XYMAX,ETAMAX,DXY,DETA;
	CHYDROtoB3D *hydrotoboltzmann;
	bool BJORKEN,COLLISIONS,INELASTIC,HYDRO_PURE_BJORKEN,DENSWRITE,BARYON_ANNIHILATION,MUTCALC;
	double ANNIHILATION_SREDUCTION;  // reduces annihilation cross section based on amount of strangeness
	int DENSWRITE_NTAU;
	int NBOSE;
	char message[200];
	double DENSWRITE_DELTAU,MUTCALC_DELTAU;
	vector<vector<vector<CMSU_BoltzmannCell *> > > cell;
	vector<vector<vector<CMuTInfo *>>> muTinfo;
	vector<double> annihilation_array;

	int Npions_fromcharges,Nprotons_fromcharges;
	
	void ReadCharges(int ichargefile);
	void GenHadronsFromCharges();
	void GenHadronsFromCharge(int balanceID,CHBCharge *charge);
	void TestChargeConservation(int pid);
	
	
	void ReadHydroInput();
	CMSUPart *GetDeadPart();
	void GetDeadParts(CMSUPart *&part1,CMSUPart *&part2);
	void GetDeadParts(array<CMSUPart*,5> &product);
	CAction *GetDeadAction();
	//vector<CLocalInfo *> localinfo;
	//
	void AnalyzeCharges(); // for testing purposes
	
	//!Constructor.
	/*!
	The constructor, which initializes all other elements of the model run.
	
	\param[in] run_name_set This is the "run name" read in from the command line.
	*/
	CMSU_Boltzmann(); // this is a constructor which does nothing but create an object
	CMSU_Boltzmann(string run_name_set,CparameterMap *parmap_set,CresList *reslist_set); // this gets all arrays ready
	CBalanceArrays *balancearrays,*bfbalancearrays;
	void CopyParMapPars(); // copies parmap to B3D variables (usually similar name)
	void InitCascade();
	void InitAnalysis();
	~CMSU_Boltzmann();	//!< Destructor.
	double tau,TAUCOLLMAX;
	int itau;
	//
	// READ IN FROM PARAMETER FILE
	long long int NACTIONSMAX,nactionstot;
	int NPARTSMAX,nbaryons,npartstot;
	int NSAMPLE,NSAMPLE_UDS2BAL;
	int DELNPARTSTOT,DELNACTIONSTOT;
	bool BINARY_RW;
	double SIGMAMAX,SIGMADEFAULT, SIGMABF,SIGMAINELASTIC, Q0; // cross sections in sq. fm
	double RESWIDTH_ALPHA; // sets spectral function for res widths
	string input_dataroot;
	string output_dataroot;
	string run_name,qualifier;
	CMSU_Decay *msudecay;

	string oscarfilename;
	FILE *oscarfile;
	//
	void SetQualifier(string qualifier_set);
	void MovePartsToFinalMap();
	double WriteOSCAR(int ievent);  // returns dnch/deta
	void ReadOSCARHeader();
	int ReadOSCAR(int ievent);
	double WriteBalanceParts(int ievent); // returns dnch/deta
	//void ReadBalanceParts();
	int ReadBalanceParts(int ievent);
	void InputPartList(CpartList *input_partlist);
	void WriteDens();
	void WriteAnnihilationData();
	
	void FindAllCollisions();
	void FindAllCellExits();
	void PerformAllActions();
	void Reset();
	void KillAllActions();
	void KillAllParts();
	void SplitPart(CMSUPart *part1,CMSUPart *part2);
	void CheckPartMap();
	void CheckActions();

	void AddAction_Activate(CMSUPart *part);
	void AddAction_Decay(CMSUPart *part,double taudecay);
	void AddAction_Collision(CMSUPart *part1,CMSUPart *part2,double tau,double pibsquared,
						double sigma_scatter,double sigma_merge,double sigma_annihilation,
						double sigma_inel,vector<double> dsigma_merge);
	void AddAction_ResetCollisions(double taureset);
	//void AddAction_SwallowParticles(double tau_breakup);
	void AddAction_ExitCell(CMSUPart *part);
	void AddAction_DensCalc(double tauwrite);
	void AddAction_MuTCalc_UpdateNPE(double taucalc);
	void AddAction_MuTCalc();

	void ListFutureCollisions();
	void PrintPartList();
	void WriteMuTInfo();
	void ReadMuTInfo();
	void WriteWeights();
	void IncrementWeightArrays();
	ChadronCount hadroncount;
	void IncrementHadronCount();
	void WriteHadronCount();
	void InitMuTCalc();

	bool FindCollision(CMSUPart *part1,CMSUPart *part2,double &taucoll);
	void Decay(CMSUPart *mother,int &nbodies,array<CMSUPart *,5> &daughter);
	void GetMassesForDecay(vector<double> &mass,int nbodies,array<CMSUPart *,5> &daughter);
	double GetSigma(CMSUPart *part1,CMSUPart *part2,double Minv2,
		double &sigma_scatter,double &sigma_merge,double &sigma_annihilation,double &sigma_inel,
		vector<double> &dsigma_merge);

	Crandy *randy;

	void PrintActionMap(CActionMap *actionmap);

	bool CheckKinematics(CMSUPart *part1,CMSUPart *part2,double &Minv2,
		double &pibsquared,double &taucoll);
	int Collide(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product,double pibsquared); // will collide if sigma>scompare

	int Collide_Scatter(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product);
	int Collide_Merge(CMSUPart *part1,CMSUPart *part2,double sigma_merge,vector<double> &dsigma_merge,
		int &nproducts,array<CMSUPart*,5> &product);
	int Collide_Annihilate(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product);
	int Collide_Inelastic(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart*,5> &product);


	void Scatter(CMSUPart *part1,CMSUPart *part2,CMSUPart *part3,CMSUPart *part4);
	bool Merge(CMSUPart *part1,CMSUPart *part2,CMSUPart *part3,CresInfo *resinfo);
	void InelasticScatter(CMSUPart *part1, CMSUPart *part2,CMSUPart *part3,CMSUPart *part4,CInelasticInfo inelinfo);
	double GetAnnihilationSigma(CMSUPart *part1,CMSUPart *part2);
	bool CancelAnnihilation(CMSUPart *part1,CMSUPart *part2);
	int Annihilate(CMSUPart *part1,CMSUPart *part2,int &nproducts,array<CMSUPart *,5> &product);

	bool ERROR_PRINT;

	long long int nscatter,nbscatter,n,nmerge,nswallow,npass,nexit;
	long long int nactivate,nannihilate,ncancel_annihilate,nregenerate,nactionkills;
	long long int ninelastic, ncollisions,oldncollisions,ndecay,ncheck,ncheck1,ncheck2;

	void freegascalc_onespecies(double m,double t,double &p,double &e,double &dens,double &sigma2,double &dedt);
	
	bool BALANCE_DECAY;
	int ibalmax;
	
	// These are used for Analysis
	CMSUPart **partarray;
	bool CALCGARRAYS;
	bool STAR_ACCEPTANCE;
	CRegenerate *regen;
	int DecayParts(int nparts);
	double CalcSpectra_PHENIX();
	void CalcSpectra_PHENIXppbar();
	double CalcSpectra_ALICE();
	void Calc3DSpectra();
	void CalcBalance();
	void CalcBalanceQGP();
	void CalcV2_STAR();
	void CalcV2_ALICE();
	void CalcHBT_STAR();
	void CalcHBT_ALICE();
	void CalcRealityDiff();
	void CalcGamma();
	void CalcMuTU();
	double legendre(int ell,double x);
	void Consolidate(string run_name);
	void CBalanceArray(CMSU_Boltzmann *boltzmann);
	void IncrementChiTotFromCharges();
	void IncrementChiTotFromHadrons();
	void DeleteCharges();
	double TotalVolume;
	Eigen::Matrix3d chitotQ,chitotH;
};

#endif
