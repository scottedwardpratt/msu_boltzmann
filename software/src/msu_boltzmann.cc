#include "boltzmann.h"
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include "boltzmann.h"
//#include "hist.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/resonances.h"
#include "msu_boltzmann/cell.h"
#include "msu_sampler/msu_sampler.h"
#include "msu_boltzmann/balancearrays.h"
#include "msu_commonutils/randy.h"
#include "msu_boltzmann/action.h"

CMSU_Boltzmann::CMSU_Boltzmann(){
};

CMSU_Boltzmann::CMSU_Boltzmann(string run_name_set){
	run_name=run_name_set;
	string parsfilename,dirname;
	dirname="model_output/"+run_name;
	parsfilename="model_output/fixed_parameters.txt";
	sprintf(message,"reading %s\n",parsfilename.c_str());
	CLog::Info(message);
	parmap.ReadParsFromFile(parsfilename);
	parsfilename=dirname+"/parameters.txt";
	sprintf(message,"reading %s\n",parsfilename.c_str());
	CLog::Info(message);
	parmap.ReadParsFromFile(parsfilename);
	CopyParMapPars();
	string command="mkdir -p model_output/"+run_name;
	system(command.c_str());
	if(BFCALC){
		parmap.ReadParsFromFile("udsdata/udsparameters.txt");
		balancearrays=new CBalanceArrays(this);
	}
	ibalmax=0;
	npartstot=nactionstot=0;
	CResList::boltzmann=this;
	CMSUPart::boltzmann=this;
	tau=0.0;
	chitotH.setZero();
	chitotQ.setZero();
	npartstot=nactionstot=0;
	PartMap.clear();
	DeadPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	randy=new CRandy(-1234);
	CAction::boltzmann=this;
	CMSUPart::boltzmann=this;
	CMSU_BoltzmannCell::boltzmann=this;
	oscarfile=NULL;
	reslist=new CResList(&parmap);

	mastersampler=new CMasterSamplre(&parmap);
	decay_nbody=new CDecay_NBody(randy);
}

void CMSU_Boltzmann::CopyParMapPars(){
	NACTIONSMAX=parmap.getI("MSU_BOLTZMANN_NACTIONSMAX",100000);
	NPARTSMAX=parmap.getI("MSU_BOLTZMANN_NPARTSMAX",200000);
	TAUCOLLMAX=parmap.getD("MSU_BOLTZMANN_TAUCOLLMAX",50.0);
	DENSWRITE=parmap.getB("MSU_BOLTZMANN_DENSWRITE",false);
	DENSWRITE_NTAU=parmap.getI("MSU_BOLTZMANN_DENSWRITE_NTAU",20);
	DENSWRITE_DELTAU=parmap.getD("MSU_BOLTZMANN_DENSWRITE_DELTAU",1.0);
	input_dataroot=parmap.getS("MSU_BOLTZMANN_INPUT_DATAROOT","udsdata/boltzmann");
	output_dataroot=parmap.getS("MSU_BOLTZMANN_OUTPUT_DATAROOT","udsdata/boltzmann");
	NSAMPLE=parmap.getI("MSU_BOLTZMANN_NSAMPLE",1);
	NSAMPLE_UDS2BAL=parmap.getI("NSAMPLE_UDS2BAL",1);
	ERROR_PRINT=parmap.getB("MSU_BOLTZMANN_ERROR_PRINT",true);
	XYMAX=parmap.getD("MSU_BOLTZMANN_XYMAX",15);
	ETAMAX=parmap.getD("MSU_BOLTZMANN_ETAMAX",1.0);
	NETA=parmap.getI("MSU_BOLTZMANN_NETA",10);
	NXY=parmap.getI("MSU_BOLTZMANN_NXY",10);
	SIGMAMAX=parmap.getD("MSU_BOLTZMANN_SIGMAMAX",30);
	NRINGSMAX=parmap.getI("MSU_BOLTZMANN_NRINGSMAX",200);
	NPRCELLSMAX=parmap.getI("MSU_BOLTZMANN_PR_NPRCELLSMAX",100000);
	COLLISIONS=parmap.getB("MSU_BOLTZMANN_COLLISIONS",true);
	BJORKEN=parmap.getB("MSU_BOLTZMANN_BJORKEN",true);
	BFCALC=parmap.getB("MSU_BOLTZMANN_BFCALC",true);
	SIGMADEFAULT=parmap.getD("MSU_BOLTZMANN_SIGMADEFAULT",1.0);
	SIGMAINELASTIC=parmap.getD( "MSU_BOLTZMANN_SIGMAINELASTIC",1.0);
	INELASTIC=parmap.getB( "MSU_BOLTZMANN_INELASTIC",false);
	NBOSE=parmap.getI("MSU_BOLTZMANN_NBOSE",1);
	Q0=parmap.getD( "MSU_BOLTZMANN_INELASTIC_Q0", 0);
	COOPERFRYE_CREATENEGPARTS=parmap.getB("MSU_BOLTZMANN_COOPERFRYE_CREATENEGPARTS","false");
	COOPERFRYE_WMAX=parmap.getI("MSU_BOLTZMANN_COOPERFRYE_WMAX",1);
	HYDRO_OCTANT_SYMMETRY=parmap.getI("HYDRO_OCTANT",2);
	HYDRO_PURE_BJORKEN=parmap.getB("HYDRO_PURE_BJORKEN",false);
	BARYON_ANNIHILATION=parmap.getB("MSU_BOLTZMANN_BARYON_ANNIHILATION",false);
	DELNPARTSTOT=parmap.getD("MSU_BOLTZMANN_DELNPARTSTOT",1000);
	DELNACTIONSTOT=parmap.getD("MSU_BOLTZMANN_DELNACTIONSTOT",2000);
	BINARY_RW=parmap.getB("MSU_BOLTZMANN_BINARY_RW",false);
	MUTCALC=parmap.getB("MSU_BOLTZMANN_MUTCALC",false);
	MUTCALC_DELTAU=parmap.getD("MSU_BOLTZMANN_MUTCALC_DELTAU",0.5);
	ANNIHILATION_SREDUCTION=parmap.getD("MSU_BOLTZMANN_ANNIHILATION_SREDUCTION",1.0);
	RESONANCE_DECAYS=parmap.getB("MSU_BOLTZMANN_RESONANCE_DECAYS",true);
	SIGMAMAX=SIGMAMAX/double(NSAMPLE);
	SIGMABF=parmap.getD("MSU_BOLTZMANN_SIGMABF",2.3);
	NPARTSMAX*=NSAMPLE;
	NACTIONSMAX*=NSAMPLE;
	DXY=XYMAX/double(NXY);
	DETA=ETAMAX/double(NETA);
	BALANCE_CALC=parmap.getB("MSU_BOLTZMANN_BALANCE_CALC",false);
}

void CMSU_Boltzmann::InitCascade(){
	// First initialize cells
	int ix,iy,ieta,jx,jy,jeta,iitau,imax;
	double xmin,xmax,ymin,ymax,etamin,etamax;
	CMSU_BoltzmannCell *c;
	CInelasticList::boltzmann=this;
	CInelasticList::UseFile = false;
	CInelasticList::UseInelasticArray = false;
	if(INELASTIC)
		inelasticlist = new CInelasticList();
	
	cell.resize(2*NXY);
	for(ix=0;ix<2*NXY;ix++){
		xmin=-XYMAX+ix*DXY;
		xmax=xmin+DXY;
		cell[ix].resize(2*NXY);
		for(iy=0;iy<2*NXY;iy++){
			ymin=-XYMAX+iy*DXY;
			ymax=ymin+DXY;
			cell[ix][iy].resize(2*NETA);
			for(ieta=0;ieta<2*NETA;ieta++){
				etamin=-ETAMAX+DETA*ieta;
				etamax=etamin+DETA;
				cell[ix][iy][ieta]=new CMSU_BoltzmannCell(xmin,xmax,ymin,ymax,etamin,etamax);
			}
		}
	}
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			for(ieta=0;ieta<2*NETA;ieta++){
				c=cell[ix][iy][ieta];
				c->ix=ix; c->iy=iy; c->ieta=ieta;
				c->creflection=NULL;
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						for(jeta=ieta-1;jeta<=ieta+1;jeta++){
							if(jx>=0 && jy>=0 && jeta>=0 && jx<2*NXY && jy<2*NXY && jeta<2*NETA){
								c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=cell[jx][jy][jeta];
							}
							else c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=NULL;
						}
					}
				}
			}
		}
	}
	if(BJORKEN){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				c=cell[ix][iy][0];
				c->ireflection=-1;
				c->creflection=cell[ix][iy][2*NETA-1];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][0]=cell[jx][jy][2*NETA-1];
						}
					}
				}				
				c=cell[ix][iy][2*NETA-1];
				c->ireflection=1;
				c->creflection=cell[ix][iy][0];
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						if(jx>=0 && jy>=0 && jx<2*NXY && jy<2*NXY){
							c->neighbor[1+jx-ix][1+jy-iy][2]=cell[jx][jy][0];
						}
					}
				}
			}
		}
	}
	if(DENSWRITE){
		for(ix=0;ix<2*NXY;ix++){
			for(iy=0;iy<2*NXY;iy++){
				for(ieta=0;ieta<2*NETA;ieta++){
					c=cell[ix][iy][ieta];
					c->dens.resize(DENSWRITE_NTAU);
					for(iitau=0;iitau<DENSWRITE_NTAU;iitau++){
						c->dens[iitau]=0.0;
					}
				}
			}
		}
	}
	if(BARYON_ANNIHILATION){
		imax=lrint(TAUCOLLMAX);
		annihilation_array.resize(imax);
		for(int i=0;i<imax;i++){
			annihilation_array[i]=0.0;
		}
	}
	if(MUTCALC || BARYON_ANNIHILATION){
		InitMuTCalc();
	}
}

void CMSU_Boltzmann::SetQualifier(string qualifier_set){
	qualifier=qualifier_set;
	string command="mkdir -p model_output/"+run_name+"/"+qualifier;
	system(command.c_str());
	if(sampler!=NULL)
		sampler->nevents=0;
	if(oscarfile!=NULL){
		fclose(oscarfile);
		oscarfile=NULL;
	}
	if(BFCALC){
		balancearrays->SetQualifier(qualifier);
	}
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.txt";
}

void CMSU_Boltzmann::Reset(){
	int iitau,ntau;
	double taucalc;
	KillAllParts();
	KillAllActions();
	tau=0.0;
	nactions=0;
	//npartstot=0;
	if(MUTCALC){
		ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
		for(iitau=1;iitau<ntau;iitau++){
			taucalc=iitau*MUTCALC_DELTAU;
			AddAction_MuTCalc_UpdateNPE(taucalc);
		}
		CMuTInfo::NETEVENTS+=NSAMPLE;
	}
}

CMSU_Boltzmann::~CMSU_Boltzmann(){
	if(oscarfile!=NULL)
		fclose(oscarfile);
}
