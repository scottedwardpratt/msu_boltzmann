#include "action.h"
#include "part.h"

void CAction::PerformCollide(){
	int colltype,iproduct,nproducts;
	double sigmatot;
	CPart *part1,*part2,*part;
	CPartMap::iterator ppos;
	CMSU_BoltzmannCell *cell;
	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	boltzmann->GetDeadParts(product);

	sigmatot=sigma_scatter+sigma_merge+sigma_annihilation+sigma_inel;
	double r=boltzmann->randy->ran();
	if(r<sigma_scatter/sigmatot){
		colltype=boltzmann->Collide_Scatter(part1,part2,nproducts,product);
	}
	else if(r<(sigma_scatter+sigma_merge)/sigmatot){
		colltype=boltzmann->Collide_Merge(part1,part2,sigma_merge,dsigma_merge,nproducts,product);
	}
	else if(r<(sigma_scatter+sigma_merge+sigma_annihilation)/sigmatot){
		colltype=boltzmann->Collide_Annihilate(part1,part2,nproducts,product);
	}
	else{
		if(!boltzmann->INELASTIC){
			sprintf(message,"In PeformCollide: Trying to Perform Inelastic Collision???");
			CLog::Fatal(message);
		}
		colltype=boltzmann->Collide_Inelastic(part1,part2,nproducts,product);
	}

	//colltype=boltzmann->Collide(part1,part2,nproducts,product,pibsquared);

	if((part1->balanceID<0 && part2->balanceID>=0) || (part1->balanceID>=0 && part2->balanceID<0)){
		colltype=-2;
	}

	if(colltype==0 || nproducts==0){
		boltzmann->npass+=1;
		part1->actionmother=boltzmann->nactions;
		part2->actionmother=boltzmann->nactions;
	}
	else if(colltype==-2){
		if(part1->balanceID>=0 && part2->balanceID<0){
			part1->CopyMomentumInfo(product[0]);
			part1->SetMass();
		}
		if(part2->balanceID>=0 && part1->balanceID<0){
			part2->CopyMomentumInfo(product[1]);
			part2->SetMass();
		}
		part1->tau_lastint=tau;
		if(part1->currentmap!=&(boltzmann->PartMap))
			part1->ChangeMap(&(boltzmann->PartMap));
		part1->active=true;
		part1->actionmother=boltzmann->nactions;
		part1->FindActions();
		part2->tau_lastint=tau;
		if(part2->currentmap!=&(boltzmann->PartMap))
			part2->ChangeMap(&(boltzmann->PartMap));
		part2->active=true;
		part2->actionmother=boltzmann->nactions;
		part2->FindActions();
	}
	else{
		part1->Kill();
		part2->Kill();
		for(iproduct=0;iproduct<nproducts;iproduct++){
			part=product[iproduct];
			part->SetMass();
			part->active=true;
			part->weight=1.0;
			part->tau_lastint=tau;
			part->actionmother=boltzmann->nactions;
			cell=part->FindCell();
			part->ChangeCell(cell);
			part->ChangeMap(&boltzmann->PartMap);
			part->FindActions();
		}
		for(iproduct=nproducts;iproduct<5;iproduct++){
			part=product[iproduct];
			part->Kill();
		}
	}
}
