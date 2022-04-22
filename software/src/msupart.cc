#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_boltzmann/msupart.h"
#include "msu_boltzmann/action.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/balancearrays.h"
#include "msu_sampler/resonances.h"
#include "msu_boltzmann/cell.h"
#include "msu_commonutils/constants.h"
#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"

CMSU_Boltzmann *CMSUPart::boltzmann=NULL;
CBalance *CMSUPart::cb=NULL;
char *CMSUPart::message=new char[500];

CMSUPart::CMSUPart(){
	cell=NULL;
	currentmap=NULL;
	bweight=1.0;
	tau0=0.0;
	currentmap=NULL;
	cell=NULL;
	active=false;
	balanceID=-999;
}

CMSUPart::CMSUPart(int keyset){
	bweight=1.0;
	key=keyset;
	tau0=0.0;
	tauexit=0.0;
	tau_lastint=0.0;
	taudecay=0.0;
	currentmap=&boltzmann->DeadPartMap;
	cell=NULL;
	actionmap.clear();
	active=false;
	resinfo=boltzmann->reslist->GetResInfoPtr(211);
	for(int alpha=0;alpha<4;alpha++){
		r[alpha]=p[alpha]=0.0;
	}
	y=eta=0.0;
	nscatt=0;
	balanceID=-999;
	boltzmann->DeadPartMap.insert(CMSUPartPair(key,this));
	boltzmann->npartstot+=1;
}

CMSUPart::~CMSUPart(){
}

void CMSUPart::Copy(CMSUPart *part){  //copies all info except actionmap
	int alpha;
	tau0=part->tau0;
	tau_lastint=part->tau_lastint;
	tauexit=part->tauexit;
	taudecay=part->taudecay;
	y=part->y;
	eta=part->eta;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
		r[alpha]=part->r[alpha];
	}
	msquared=part->msquared;
	resinfo=part->resinfo;
	balanceID=part->balanceID;
	bweight=part->bweight;
	//eta0=part->eta0;
	phi0=part->phi0;
}

void CMSUPart::CopyPositionInfo(CMSUPart *part){  //copies all info except actionmap
	int alpha;
	tau0=part->tau0;
	eta=part->eta;
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=part->r[alpha];
	}
}

void CMSUPart::CopyMomentumInfo(CMSUPart *part){  //copies all info except actionmap
	int alpha;
	y=part->y;
	msquared=part->msquared;
	for(alpha=0;alpha<4;alpha++){
		p[alpha]=part->p[alpha];
	}
	bweight=part->bweight;
}

void CMSUPart::InitBalance(int IDset,double rxset,double ryset,double tauset,double etaset,double pxset,double pyset,double mset,double rapidityset,double bweightset,int balanceIDset){
	Init(IDset,rxset,ryset,tauset,etaset,pxset,pyset,mset,rapidityset,bweightset);
	balanceID=balanceIDset;
}

void CMSUPart::Init(int IDset,double rxset,double ryset,double tauset,double etaset,double pxset,double pyset,double mset,double rapidityset,double bweightset){
	double et;
	CresInfo *resinfoptr;
	int ID;
	resinfo=boltzmann->reslist->GetResInfoPtr(IDset);
	ID=resinfo->pid;
	if(ID!=IDset){
		sprintf(message,"ID mismatch, ID=%d, resinfo->pidID=%d\n",IDset,ID);
		CLog::Info(message);
	}
	p[1]=pxset; p[2]=pyset; msquared=mset*mset; y=rapidityset;
	r[1]=rxset; r[2]=ryset; tau0=tauset; eta=etaset;
	//eta0=eta;
	phi0=atan2(r[2],r[1]);
	tau_lastint=tau0-1.0E-6;
	nscatt=0;
	bweight=bweightset;
	if(fabs(eta)>boltzmann->ETAMAX && boltzmann->BJORKEN){
		sprintf(message,"in part->init, eta out of bounds, =%g\n",eta);
		CLog::Info(message);
		sprintf(message,"boltzmann->ETAMAX=%g\n",boltzmann->ETAMAX);
		CLog::Fatal(message);
	}
	resinfoptr=boltzmann->reslist->GetResInfoPtr(ID);
	if(resinfoptr->decay==false){
		msquared=resinfoptr->mass*resinfoptr->mass;
	}
	r[3]=tau0*sinh(eta);
	r[0]=tau0*cosh(eta);
	et=sqrt(p[1]*p[1]+p[2]*p[2]+msquared);
	p[3]=et*sinh(y);
	Setp0();
	if(tau0<0.0){
		Print();
		sprintf(message,"FATAL: tau0<0, tau0^2=%g\n",tau0);
		CLog::Fatal(message);
	}
	if(boltzmann->BJORKEN && fabs(eta)>boltzmann->ETAMAX){
		CyclicReset();
		sprintf(message,"performed cyclic reset in CMSUPart::Init()\n");
		CLog::Fatal(message);
	}
	active=false;
	ChangeMap(&(boltzmann->PartMap));
	boltzmann->AddAction_Activate(this);
	actionmother=boltzmann->nactions;
}

void CMSUPart::CheckRapidity(){
	if(fabs(y-atanh(p[3]/p[0]))>0.001){
		Print();
		sprintf(message,"rapidity screwed up!\n");
		CLog::Fatal(message);
	}
}

void CMSUPart::CyclicReset(){
	double eta_offset,etamax=boltzmann->ETAMAX;
	double mt;
	while(fabs(eta)>etamax){
		eta_offset=2.0*etamax;
		if(eta<0.0)
			eta_offset=-eta_offset;
		eta-=eta_offset;
		y-=eta_offset;
		mt=sqrt(msquared+p[1]*p[1]+p[2]*p[2]);
		p[3]=mt*sinh(double(y));
		p[0]=sqrt(mt*mt+p[3]*p[3]);
		r[3]=tau0*sinh(double(eta));
		r[0]=tau0*cosh(double(eta));
	}
}

void CMSUPart::Print(){
	sprintf(message,"________________ PART INFO FOR key=%d _____________________________\n",key);
	CLog::Info(message);
	sprintf(message,"Minv^2=%g, ID=%d -----  %s ------\n",p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3],resinfo->pid,resinfo->name.c_str());
	CLog::Info(message);
	sprintf(message,"ID=%d, m_onshell=%g, M=%g, tau0=%g=?%g, tauexit=%g\n r=(%g,%g,%g,%g) eta=%g=?%g\n", 
	resinfo->pid,resinfo->mass,sqrt(msquared),double(tau0),sqrt(r[0]*r[0]-r[3]*r[3]),tauexit,r[0],r[1],r[2],r[3],eta,GetEta(tau0));
	CLog::Info(message);
	sprintf(message,"bweight=%g, key=%d, actionmother=%d, active=%d, balanceID=%d\n",
	bweight,key,actionmother,int(active),balanceID);
	CLog::Info(message);
	sprintf(message,"p=(%15.9e,%15.9e,%15.9e,%15.9e), y=%g =? %g\n",p[0],p[1],p[2],p[3],y,atanh(p[3]/p[0]));
	CLog::Info(message);
	string currentmapname="IN CELL";
	if(currentmap==&(boltzmann->PartMap)) currentmapname="PartMap";
	if(currentmap==&(boltzmann->DeadPartMap)) currentmapname="DeadPartMap";
	sprintf(message,"currentmap=%s\n",currentmapname.c_str());
	CLog::Info(message);
	if(cell==NULL) sprintf(message,"CELL=NULL\n");
	CLog::Info(message);
	if(nextcell==NULL) sprintf(message,"NEXTCELL=NULL\n");
	CLog::Info(message);
	if(cell!=NULL) sprintf(message,"Cell No: ix=%d, iy=%d, ieta=%d\n",cell->ix,cell->iy,cell->ieta);
	CLog::Info(message);
	sprintf(message,"________________________________________________________________________\n");
	CLog::Info(message);CLog::Info(message);
}

void CMSUPart::ChangeMap(CMSUPartMap *newmap){
	if(newmap!=currentmap){
		DeleteFromCurrentMap();
		AddToMap(newmap);
	}
}

CMSUPartMap::iterator CMSUPart::DeleteFromCurrentMap(){
	CMSUPartMap::iterator neighbor;
	CMSUPartMap::iterator ppos=GetPos(currentmap);
	neighbor=ppos;
	neighbor++;
	if(ppos==currentmap->end()){
		sprintf(message,"FATAL: In CMSUPart::DeleteFromCurrentMap, can't find ppos!!!\n");
		Print();
		sprintf(message,"currentmap has length %d\n",int(currentmap->size()));
		CLog::Fatal(message);
	}
	else currentmap->erase(ppos);
	currentmap=NULL;
	return neighbor;
}

CMSUPartMap::iterator CMSUPart::DeleteFromMap(CMSUPartMap *partmap){
	CMSUPartMap::iterator neighbor;
	CMSUPartMap::iterator ppos=GetPos(partmap);
	neighbor=ppos;
	if(ppos==partmap->end()){
		Print();
		sprintf(message,"FATAL: In CMSUPart::DeleteFromMap, can't find ppos!!!\n");
		CLog::Fatal(message);
	}
	else{
		neighbor++;
		partmap->erase(ppos);
	}
	//partmap=NULL;
	return neighbor;
}

void CMSUPart::AddToMap(CMSUPartMap *newmap){
	newmap->insert(CMSUPartPair(key,this));
	currentmap=newmap;
}

void CMSUPart::AddToMap(CMSUPartMap::iterator guess,CMSUPartMap *newmap){
	newmap->insert(guess,CMSUPartPair(key,this));
	currentmap=newmap;
}

void CMSUPart::SubtractAction(CAction *action){
	CActionMap::iterator epos=action->GetPos(&actionmap);
	if(epos!=actionmap.end())
		actionmap.erase(epos);
}

void CMSUPart::AddAction(CAction *action){
	actionmap.insert(CActionPair(action->key,action));
}

void CMSUPart::Propagate(double tau){
	if(boltzmann->BJORKEN && fabs(eta)>boltzmann->ETAMAX){
		sprintf(message,"eta screwy before propagation\n");
		CLog::Info(message);
		sprintf(message,"eta=%g\n",eta);
		CLog::Info(message);
	}
	double t0;
	CMSUPartMap::iterator neighbor;
	if(active==true){
		eta=GetEta(tau);//y-asinh((tau0/tau)*sinh(y-eta));
		if(currentmap==&(boltzmann->PartMap) && boltzmann->tau<boltzmann->TAUCOLLMAX && fabs(eta)>boltzmann->ETAMAX && boltzmann->BJORKEN && boltzmann->COLLISIONS){
			//sprintf(message,"eta out of bounds after propagation,correcting, etai=%g, etaf=%g, taui=%g, tauf=%g\n",etai,eta,tau0,tau);
			//CLog::Info(message);
			//Print();
			if(eta>boltzmann->ETAMAX)
				eta-=2.0*boltzmann->ETAMAX;
			if(eta<-boltzmann->ETAMAX)
				eta+=2.0*boltzmann->ETAMAX;
			r[0]=tau0*cosh(eta);
			r[3]=tau0*sinh(eta);
			//Misc::Pause();
		}
		tau0=tau;
		t0=r[0];
		r[0]=tau0*cosh(eta);
		r[3]=tau0*sinh(eta);
		r[1]+=(p[1]/p[0])*(r[0]-t0);
		r[2]+=(p[2]/p[0])*(r[0]-t0);
	}
	else{
		r[0]=tau*cosh(eta);
		r[3]=tau*sinh(eta);
		tau0=tau;
	}
}

CMSUPartMap::iterator CMSUPart::GetPos(CMSUPartMap *pmap){
	CMSUPartMap::iterator ppos=pmap->find(key);
	return ppos;
}

void CMSUPart::CheckMap(CMSUPartMap *expectedpartmap){
	if(currentmap!=expectedpartmap){
		
		Print();
		if(currentmap==&(boltzmann->DeadPartMap)){
			sprintf(message,"particle in DeadPartMap\n");
			CLog::Info(message);
		}
		if(currentmap==&(boltzmann->PartMap)){
			sprintf(message,"particlein PartMap\n");
			CLog::Info(message);
		}
		sprintf(message,"XXXXXXXXX particle not in expected map XXXXXXXXX\n");
		CLog::Fatal(message);
	}
}

double CMSUPart::GetMass(){
	if(resinfo->pid==22)
		return 0.0;
	else
		return sqrt(msquared);
}

/* in sampler there is an array of densities and temperature. also make array of densities of minmass*/

void CMSUPart::Setp0(){
	p[0]=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]+msquared);
}

void CMSUPart::SetY(){
	y=asinh(p[3]/GetMT());
}

void CMSUPart::SetMass(){
	msquared=p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3];
}

double CMSUPart::GetEta(double tau){
	double dy,deta,dtau0,dtau;
	if(active){
		dy=y;
		deta=eta;
		dtau0=tau0;
		dtau=tau;
		deta=dy-asinh((dtau0/dtau)*sinh(dy-deta));
		return deta;
	}
	else return eta;
}

double CMSUPart::GetMT(){
	if(p[0]<fabs(p[3])){
		Print();
		sprintf(message,"CMSUPart::GetMT, catastrophe\n");
		CLog::Fatal(message);
	}
	return sqrt(p[0]*p[0]-p[3]*p[3]);
}

void CMSUPart::KillActions(){
	CActionMap::iterator ep,epp;
	CAction *action;
	ep=actionmap.begin();
	while(ep!=actionmap.end()){
		action=ep->second;
		epp=ep; ++epp;
		action->Kill();
		ep=epp;
	}
	actionmap.clear();
}

void CMSUPart::Kill(){
	KillActions();
	if(cell!=NULL){
		RemoveFromCell();
		cell=NULL;
	}
	DeleteFromCurrentMap();
	eta=0.0;
	nscatt=0;
	tau0=tau_lastint=tauexit=-1.0;
	AddToMap(boltzmann->DeadPartMap.begin(),&(boltzmann->DeadPartMap));
	bweight=1.0;
	active=false;
	balanceID=-1;
}

void CMSUPart::BjorkenTranslate(){
	if(eta<-boltzmann->ETAMAX || eta>boltzmann->ETAMAX){
		sprintf(message,"eta out of bounds before translation\n");
		CLog::Info(message);
		Print();
		if(eta>boltzmann->ETAMAX)
			eta=boltzmann->ETAMAX;
		if(eta<-boltzmann->ETAMAX)
			eta=-boltzmann->ETAMAX;
		r[0]=tau0*cosh(eta);
		r[3]=tau0*sinh(eta);
	}
	double mt;
	if(cell->ieta==0){
		eta+=2.0*boltzmann->ETAMAX;
		y+=2.0*boltzmann->ETAMAX;
	}
	else{
		eta-=2.0*boltzmann->ETAMAX;
		y-=2.0*boltzmann->ETAMAX;
	}
	mt=sqrt(p[0]*p[0]-p[3]*p[3]);
	p[3]=mt*sinh(y);
	p[0]=mt*cosh(y);
	r[0]=tau0*cosh(eta);
	r[3]=tau0*sinh(eta);
	SetMass();
}

void CMSUPart::BjorkenUnTranslate(){
	double mt;
	if(eta>boltzmann->ETAMAX-1.0E-10){
		eta-=2.0*boltzmann->ETAMAX;
		y-=2.0*boltzmann->ETAMAX;
	}
	else{
		eta+=2.0*boltzmann->ETAMAX;
		y+=2.0*boltzmann->ETAMAX;
	}
	mt=sqrt(p[0]*p[0]-p[3]*p[3]);
	p[3]=mt*sinh(y);
	p[0]=mt*cosh(y);
	r[0]=tau0*cosh(eta);
	r[3]=tau0*sinh(eta);
	SetMass();
}

void CMSUPart::FindCollisions(){
	int ix,iy,ieta;
	double taucoll;
	CMSUPart *part2,*part1=this;
	CMSUPartMap::iterator ppos;
	CMSU_BoltzmannCell *cell2;

	for(ix=0;ix<3;ix++){
		for(iy=0;iy<3;iy++){
			for(ieta=0;ieta<3;ieta++){
				cell2=cell->neighbor[ix][iy][ieta];
				if(cell2!=NULL){
					ppos=cell2->partmap.begin();
					while(ppos!=cell2->partmap.end()){
						part2=ppos->second;
						if(part1!=part2 && part1->actionmother!=part2->actionmother){
							if(part1->balanceID<0 || part2->balanceID<0){
								boltzmann->FindCollision(part1,part2,taucoll);
							}
						}
						++ppos;
					}
				}
			}
		}
	}
}

CMSU_BoltzmannCell *CMSUPart::FindCell(){
	if(tau0>boltzmann->TAUCOLLMAX || !boltzmann->COLLISIONS){
		return NULL;
	}
	int ieta,ix,iy;
	double deta=boltzmann->ETAMAX/double(boltzmann->NETA);
	ieta=lrint(floor((eta+boltzmann->ETAMAX)/deta));
	if(ieta<0 ||ieta>=2*boltzmann->NETA){
		return NULL;
	}
	double dx=boltzmann->XYMAX/double(boltzmann->NXY);
	double dy=dx;
	ix=lrint(floor((r[1]+boltzmann->XYMAX)/dx));
	if(ix<0 || ix>=2*boltzmann->NXY){
		return NULL;
	}
	iy=lrint(floor((r[2]+boltzmann->XYMAX)/dy));
	if(iy<0 || iy>=2*boltzmann->NXY){
		return NULL;
	}
	return boltzmann->cell[ix][iy][ieta];
}

void CMSUPart::FindDecay(){
	double t,gamma,vz,newt,newz;
	t=HBARC_GEV/resinfo->width;
	gamma=p[0]/sqrt(msquared);
	t=-t*gamma*log(boltzmann->randy->ran());
	vz=p[3]/p[0];
	newt=r[0]+t;
	newz=r[3]+vz*t;
	taudecay=sqrt(newt*newt-newz*newz);
	if(taudecay<tauexit || tauexit>boltzmann->TAUCOLLMAX || cell==NULL){
		boltzmann->AddAction_Decay(this,taudecay);
	}
}

void CMSUPart::FindCellExit(){
	if(active){
		double t,taux,tauy,taueta,z;
		double etamax=cell->etamax,etamin=cell->etamin;
		nextcell=NULL;
		tauexit=1.0E50;
		taueta=taux=tauy=tauexit;
		double vx=p[1]/p[0];
		double vy=p[2]/p[0];
		double vz=p[3]/p[0];

		if(vx>0)
			t=(cell->xmax-r[1])/vx;
		else
			t=(cell->xmin-r[1])/vx;
		t=t+r[0];
		z=r[3]+vz*(t-r[0]);
		taux=sqrt(t*t-z*z);

		if(vy>0)
			t=(cell->ymax-r[2])/vy;
		else
			t=(cell->ymin-r[2])/vy;
		t=t+r[0];
		z=r[3]+vz*(t-r[0]);
		tauy=sqrt(t*t-z*z);

		if(y<cell->etamin)
			taueta=tau0*sinh(y-eta)/sinh(y-etamin);
		else if(y>cell->etamax)
			taueta=tau0*sinh(y-eta)/sinh(y-etamax);

		if(taux<tauy && taux<taueta){
			tauexit=taux;
			if(vx<0)
				nextcell=cell->neighbor[0][1][1];
			else
				nextcell=cell->neighbor[2][1][1];
		}
		else if (tauy<taueta){
			tauexit=tauy;
			if(vy<0)
				nextcell=cell->neighbor[1][0][1];
			else
				nextcell=cell->neighbor[1][2][1];
		}
		else{
			tauexit=taueta;
			if(y<etamin)
				nextcell=cell->neighbor[1][1][0];
			else
				nextcell=cell->neighbor[1][1][2];
		}
		if(tauexit<boltzmann->TAUCOLLMAX)
			boltzmann->AddAction_ExitCell(this);
	}
}

void CMSUPart::FindActions(){
	KillActions();
	if(active!=true){
		sprintf(message,"CMSUPart::FindActions(), trying to Reset Inactive particle\n");
		CLog::Info(message);
		KillActions();
	}
	if(resinfo->pid!=22 && msquared<resinfo->minmass*resinfo->minmass-1.0E-6){
		sprintf(message,"msquared too small, M=%14.9e, minmass=%14.9e\n",sqrt(msquared),resinfo->minmass);
		CLog::Info(message);
		Print();
		exit(1);
		msquared=resinfo->minmass*resinfo->minmass;
		Setp0();
	}

	if(cell!=NULL){
		if(boltzmann->COLLISIONS && tau0<boltzmann->TAUCOLLMAX){
			FindCellExit();
			FindCollisions();
		}
		else{
			ChangeCell(NULL);
		}
	}
	if(resinfo->decay)
		FindDecay();
	if(currentmap!=&(boltzmann->PartMap)){
		sprintf(message,"FindActions: part in wrong map\n");
		Print();
		CLog::Fatal(message);
	}
}

double CMSUPart::GetPseudoRapidity(){
	double pmag,eta_ps;
	pmag=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
	eta_ps=atanh(p[3]/pmag);
	return eta_ps;
}

double CMSUPart::GetRapidity(){
	return 0.5*log((p[0]+p[3])/(p[0]-p[3]));
}

void CMSUPart::Boost(FourVector &u){
	BoostP(u);
	BoostR(u);
}

void CMSUPart::BoostP(FourVector &u){
	int alpha;
	FourVector pprime;
	Misc::Boost(u,p,pprime);
	for(alpha=0;alpha<4;alpha++)
		p[alpha]=pprime[alpha];
	y=atanh(p[3]/p[0]);
}

void CMSUPart::BoostR(FourVector &u){
	int alpha;
	FourVector rprime;
	Misc::Boost(u,r,rprime);
	for(alpha=0;alpha<4;alpha++)
		r[alpha]=rprime[alpha];
	eta=atanh(r[3]/r[0]);
	tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
}

void CMSUPart::RemoveFromCell(){
	if(cell!=NULL){
		CMSUPartMap::iterator ppos=GetPos(&(cell->partmap));
		if(ppos==cell->partmap.end()){
			sprintf(message,"FATAL: In CMSUPart::RemoveFromCell, can't find ppos!!!\n");
			Print();
			sprintf(message,"cell partmap has length %d\n",int(cell->partmap.size()));
			CLog::Fatal(message);
		}
		else{
			cell->partmap.erase(ppos);
		}
	}
}

void CMSUPart::CheckCell(){
	if(cell!=NULL){
		CMSUPartMap::iterator ppos=GetPos(&(cell->partmap));
		if(ppos==cell->partmap.end()){
			sprintf(message,"FATAL: In CMSUPart::RemoveFromCell, can't find ppos!!!\n");
			Print();
			sprintf(message,"cell partmap has length %d\n",int(cell->partmap.size()));
			CLog::Fatal(message);
		}
	}
}

void CMSUPart::ChangeCell(CMSU_BoltzmannCell *newcell){
	if(newcell!=cell){
		if(cell!=NULL)
			RemoveFromCell();
		if(newcell!=NULL){
			newcell->partmap.insert(CMSUPartPair(key,this));
		}
		cell=newcell;
	}
}

void CMSUPart::GetHBTPars(double &t,double &rout,double &rside,double &rlong){
	const double tcompare=15.0;
	double pt,ptsquared,et;
	rlong=tau0*sinh(eta-y);
	t=tau0*cosh(eta-y);
	ptsquared=p[1]*p[1]+p[2]*p[2];
	pt=sqrt(ptsquared);
	et=sqrt(ptsquared+msquared);
	rout=(p[1]*r[1]+p[2]*r[2])/pt;
	rout=rout-(pt/et)*(t-tcompare);
	rside=(p[1]*r[2]-p[2]*r[1])/pt;
}

void CMSUPart::BoostRap(double dely){
	double gamma=cosh(dely),gammav=sinh(dely);
	double r0=r[0],p0=p[0];
	y=y+dely;
	eta=eta+dely;
	r[0]=gamma*r0+gammav*r[3];
	r[3]=gamma*r[3]+gammav*r0;
	p[0]=gamma*p0+gammav*p[3];
	p[3]=gamma*p[3]+gammav*p0;
}

void CMSUPart::CalcDCA(double *dca){
	int alpha;
	char nantestc[20];
	string nantests;
	double pdotr,p2;
	p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
	pdotr=(p[1]*r[1]+p[2]*r[2]+p[3]*r[3])/p2;
	for(alpha=1;alpha<4;alpha++)
		dca[alpha]=(r[alpha]-pdotr*p[alpha])/1.0E13;
	dca[0]=sqrt(dca[1]*dca[1]+dca[2]*dca[2]+dca[3]*dca[3]);
	sprintf(nantestc,"%g",r[0]);
	nantests=nantestc;
	if(nantests=="NaN" || nantests=="nan" || nantests=="inf" || nantests=="INF"){
		sprintf(message,"::: dca=(%g,%g,%g,%g)\n",dca[0],dca[1],dca[2],dca[3]);
		sprintf(message,"::: r=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
		sprintf(message,"::: p=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
		CLog::Fatal(message);
	}
}