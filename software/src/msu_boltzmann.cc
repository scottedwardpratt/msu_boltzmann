#include "msu_boltzmann/msu_boltzmann.h"
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <filesystem>
#include "msu_boltzmann/msupart.h"
#include "msu_eos/resonances.h"
#include "msu_boltzmann/cell.h"
#include "msu_sampler/sampler.h"
#include "msu_boltzmann/balancearrays.h"
#include "msu_commonutils/randy.h"
#include "msu_boltzmann/action.h"

CmasterSampler *CMSU_Boltzmann::mastersampler=nullptr;

using namespace std;
using namespace NMSUPratt;

CMSU_Boltzmann::CMSU_Boltzmann(){
};

CMSU_Boltzmann::CMSU_Boltzmann(int run_number,int subrun_number_set,CresList *reslist_set){
	run_name="run"+to_string(run_number);
	subrun_number=subrun_number_set;
	reslist=reslist_set;
	string parsfilename="modelruns/fixed_parameters.txt";
	parmap.ReadParsFromFile(parsfilename);
	parsfilename="modelruns/run"+to_string(run_number)+"/mod_parameters.txt";
	if(std::filesystem::exists(parsfilename))
		parmap.ReadParsFromFile(parsfilename);
	nevents=0;
	CopyParMapPars();
	if(BFCALC){
		balancearrays=new CBalanceArrays(this);
	}
	ibalmax=0;
	npartstot=nactionstot=0;
	//CresList::boltzmann=this;
	CMSUPart::boltzmann=this;
	tau=0.0;
	chargevec.resize(40000);
	allpartsvec.resize(DELNPARTSTOT);
	allactionsvec.resize(DELNACTIONSTOT);
	chitotH.setZero();
	chitotQ.setZero();
	PartMap.clear();
	DeadPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	randy=new Crandy(run_number*1000+subrun_number);
	msudecay=new CMSU_Decay(randy);
	CAction::boltzmann=this;
	for(int iaction=0;iaction<DELNACTIONSTOT;iaction++){
		allactionsvec[iaction]=new CAction(iaction);
	}
	for(int ipart=0;ipart<DELNPARTSTOT;ipart++){
		allpartsvec[ipart]=new CMSUPart(ipart);
	}
	CMSUPart::boltzmann=this;
	CMSU_BoltzmannCell::boltzmann=this;
	msudecay->boltzmann=this;
	oscarfile=nullptr;
}

void CMSU_Boltzmann::CopyParMapPars(){
	NPARTSMAX=parmap.getI("MSU_BOLTZMANN_NPARTSMAX",200000);
	TAUCOLLMAX=parmap.getD("MSU_BOLTZMANN_TAUCOLLMAX",50.0);
	DENSWRITE=parmap.getB("MSU_BOLTZMANN_DENSWRITE",false);
	DENSWRITE_NTAU=parmap.getI("MSU_BOLTZMANN_DENSWRITE_NTAU",20);
	DENSWRITE_DELTAU=parmap.getD("MSU_BOLTZMANN_DENSWRITE_DELTAU",1.0);
	input_dataroot=parmap.getS("MSU_BOLTZMANN_INPUT_DATAROOT","modelruns");
	output_dataroot=parmap.getS("MSU_BOLTZMANN_OUTPUT_DATAROOT","modelruns");
	NSAMPLE=parmap.getI("MSU_SAMPLER_NSAMPLE",1);
	NSAMPLE_UDS2BAL=parmap.getI("MSU_BOLTZMANN_NSAMPLE_UDS2BAL",1);
	ERROR_PRINT=parmap.getB("MSU_BOLTZMANN_ERROR_PRINT",true);
	XYMAX=parmap.getD("MSU_BOLTZMANN_XYMAX",15);
	ETAMAX=parmap.getD("MSU_BOLTZMANN_ETAMAX",1.0);
	parmap.set("MSU_SAMPLER_BJORKEN_ETAMAX",ETAMAX);
	NETA=parmap.getI("MSU_BOLTZMANN_NETA",10);
	NXY=parmap.getI("MSU_BOLTZMANN_NXY",10);
	COLLISIONS=parmap.getB("MSU_BOLTZMANN_COLLISIONS",true);
	BJORKEN=parmap.getB("MSU_BOLTZMANN_BJORKEN",true);
	BFCALC=parmap.getB("MSU_BOLTZMANN_BFCALC",false);
	SIGMADEFAULT=parmap.getD("MSU_BOLTZMANN_SIGMADEFAULT",1.0);
	SIGMAINELASTIC=parmap.getD("MSU_BOLTZMANN_SIGMAINELASTIC",1.0);
	SIGMABF=parmap.getD("MSU_BOLTZMANN_SIGMABF",2.3);
	INELASTIC=parmap.getB( "MSU_BOLTZMANN_INELASTIC",false);
	NBOSE=parmap.getI("MSU_BOLTZMANN_NBOSE",1);
	Q0=parmap.getD( "MSU_BOLTZMANN_INELASTIC_Q0", 0);
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
	DELPT_SPECTRA=parmap.getD("MSU_BOLTZMANN_DELPT_SPECTRA",0.025);
	DELPT_V2=parmap.getD("MSU_BOLTZMANN_DELPT_V2",0.025);
	NPT_SPECTRA=parmap.getI("MSU_BOLTZMANN_NPT_SPECTRA",120);
	NPT_V2=parmap.getI("MSU_BOLTZMANN_NPT_V2",120);
	CHARGEVEC_MAXSIZE=parmap.getI("MSU_BOLTZMANN_CHARGEMAP_MAXSIZE",40000);
	
	NPARTSMAX*=NSAMPLE;
	DXY=XYMAX/double(NXY);
	DETA=ETAMAX/double(NETA);
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
				c->creflection=nullptr;
				for(jx=ix-1;jx<=ix+1;jx++){
					for(jy=iy-1;jy<=iy+1;jy++){
						for(jeta=ieta-1;jeta<=ieta+1;jeta++){
							if(jx>=0 && jy>=0 && jeta>=0 && jx<2*NXY && jy<2*NXY && jeta<2*NETA){
								c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=cell[jx][jy][jeta];
							}
							else c->neighbor[1+jx-ix][1+jy-iy][1+jeta-ieta]=nullptr;
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
	meanpt.resize(3);
	meanv2.resize(3);
	meanpt_denom.resize(3);
	meanv2_denom.resize(3);
	spectra.resize(3);
	for(int ispecies=0;ispecies<3;ispecies++){
		meanpt[ispecies]=0.0;
		meanv2[ispecies]=0.0;
		meanpt_denom[ispecies]=0.0;
		meanv2_denom[ispecies]=0.0;
		spectra[ispecies].resize(NPT_SPECTRA);
		for(int ipt=0;ipt<NPT_V2;ipt++)
			spectra[ispecies][ipt]=0.0;
	}
	v2.resize(3);
	v2denom.resize(3);
	for(int ispecies=0;ispecies<3;ispecies++){
		v2[ispecies].resize(NPT_V2);
		v2denom[ispecies].resize(NPT_V2);
		for(int ipt=0;ipt<NPT_V2;ipt++){
			v2[ispecies][ipt]=0.0;
			v2denom[ispecies][ipt]=0.0;
		}
	}
	if(MUTCALC || BARYON_ANNIHILATION){
		InitMuTCalc();
	}
}

void CMSU_Boltzmann::SetQualifier(string qualifier_set){
	qualifier=qualifier_set;
	string command="mkdir -p modelruns/"+run_name+"/"+qualifier;
	system(command.c_str());
	if(oscarfile!=nullptr){
		fclose(oscarfile);
		oscarfile=nullptr;
	}
	if(BFCALC){
		balancearrays->SetQualifier(qualifier);
	}
	oscarfilename="modelruns/"+run_name+"/"+qualifier+"/oscar.txt";
}

void CMSU_Boltzmann::Reset(){
	int iitau,ntau;
	double taucalc;
	KillAllActions();
	KillAllParts();
	PartMap.clear();
	DeadPartMap.clear();
	ActionMap.clear();
	DeadActionMap.clear();
	for(unsigned int ipart=0;ipart<allpartsvec.size();ipart++){
		allpartsvec[ipart]->InitDead(int(ipart));
	}
	for(unsigned int iaction=0;iaction<allactionsvec.size();iaction++){
		allactionsvec[iaction]->InitDead(int(iaction));
	}
	tau=0.0;
	nactionstot=0;
	npartstot=0;
	if(MUTCALC){
		ntau=lrint(TAUCOLLMAX/MUTCALC_DELTAU);
		for(iitau=0;iitau<ntau;iitau++){
			taucalc=(iitau+1)*MUTCALC_DELTAU;
			AddAction_MuTCalc_UpdateNPE(taucalc);
		}
		CMuTInfo::NETEVENTS+=NSAMPLE;
	}
	randy->reset(nevents);
	nevents+=1;
}

CMSU_Boltzmann::~CMSU_Boltzmann(){
	if(oscarfile!=nullptr)
		fclose(oscarfile);
}
