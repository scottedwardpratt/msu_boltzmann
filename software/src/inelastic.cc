#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_boltzmann/inelastic.h"
#include "msu_boltzmann/boltzmanndefs.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_eos/resonances.h"
#include "msu_commonutils/misc.h"

using namespace std;
using namespace NMSUPratt;

CMSU_Boltzmann *CInelasticList::boltzmann=NULL;
bool CInelasticList::UseFile = false;
bool CInelasticList::UseInelasticArray = false;

CInelasticList::CInelasticList(){
	ifstream inelasticfile;
	// UseFile = false;
	// UseInelasticArray = true;

	if(boltzmann!=NULL){
		NResonances = boltzmann->reslist->resmap.size();
		filename = boltzmann->parmap.getS("MSU_BOLTZMANN_INELASTIC_INFO_FILE",string("inelastic.tmp"));

		//create a new NResonances x NResonances array to store lists of elements
		InelasticArray = new list<CInelasticInfo> *[NResonances];
		for(int i = 0; i<NResonances; i++){
			InelasticArray[i] = new list<CInelasticInfo> [NResonances];
		}
		if(UseFile){
			inelasticfile.open(filename.c_str());
			if(inelasticfile){
				ReadInelasticInfo(true);
			}
			else{
				ReadInelasticInfo(false);
			}
		}
		else{
			ReadInelasticInfo(false);
		}
	}
}

CInelasticInfo::CInelasticInfo(CresInfo *resinfo_1_in, CresInfo *resinfo_2_in, int type_in){
	//by convention, resinfo_1 will have the lighter of the two particles
	if(resinfo_1_in->mass <= resinfo_2_in->mass){
		resinfo_1 = resinfo_1_in;
		resinfo_2 = resinfo_2_in;
	}
	else{
		resinfo_1 = resinfo_2_in;
		resinfo_2 = resinfo_1_in;
	}
	type = type_in;
	min_mass = resinfo_1->mass + resinfo_2->mass;
	net_q = resinfo_1->charge + resinfo_2->charge;
	net_s = resinfo_1->strange + resinfo_2->strange;
	net_b = resinfo_1->baryon + resinfo_2->baryon;
}

void CInelasticInfo::Print(){
	resinfo_1->Print();
	resinfo_2->Print();
}

void CInelasticList::ReadInelasticInfo(bool FromFile){
	/*
	NOTE: The thermal array indices are named the convention {absolute value of quantities}{sign of quantites}
	and the quantities are named in alphabetical order; baryon number (b), charge (q), strangeness (s).
	The sign function returns 0 if the number is negative, and 1 otherwise.
	*/

	string filler;
	fstream inelasticfile;
	int ires1, ires2, ires3, ires4, netq, netb, nets, pmq=0, pmb=0, pms=0, size, sum = 0;
	double foo = 0;
	CresInfo *resinfoptr_1 = NULL,*resinfoptr_2 = NULL;
	CresInfoMap::iterator rpos1,rpos2;
	CInelasticInfo *temp;
	list<CInelasticInfo>::iterator Th_iter;
	CresList *reslist=boltzmann->reslist;
	int nres=reslist->resmap.size();

	if(FromFile){
		inelasticfile.open(filename.c_str(), fstream::in);
		if(inelasticfile){
			inelasticfile >> filler; //catches #THERMAL tag
			while(sum != 18 && inelasticfile.good() ){
				inelasticfile >> netb >> netq >> nets >> pmb >> pmq >> pms >> size;
				sum = netq + netb + nets + pmq + pmb + pms;
				for(int i=0; i<size; i++){
					inelasticfile >> ires1 >> ires2;
					resinfoptr_1=boltzmann->reslist->GetResInfoPtr(ires1);
					resinfoptr_2=boltzmann->reslist->GetResInfoPtr(ires2);

					CInelasticInfo temp2(resinfoptr_1,resinfoptr_2, 0);
					ThermalArray[netb][netq][nets][pmb][pmq][pms].push_back(temp2);
				}
			}
			if(UseInelasticArray){
				inelasticfile >> filler;
				if(strcmp(filler.c_str(), "#INELASTIC_ARRAY") == 0){
					while(inelasticfile >> ires1 >> ires2 >> size){
						for(int i=0; i<size; i++){
							inelasticfile >>ires3 >> ires4;
							resinfoptr_1=boltzmann->reslist->GetResInfoPtr(ires3);
							resinfoptr_2=boltzmann->reslist->GetResInfoPtr(ires4);

							CInelasticInfo temp2(resinfoptr_1,resinfoptr_2,0);
							InelasticArray[ires1][ires2].push_back(temp2);
						}
					}
				}else{
					inelasticfile.close();
					remove(filename.c_str());
					ReadInelasticInfo(false);
					return;
				}
			}
		}
		else {
			CLog::Fatal("Error; inelastic file can't be opened\n");
		}
		inelasticfile.close();
	}
	else{
		for(rpos1=reslist->resmap.begin();rpos1!=reslist->resmap.end();rpos1++){
			resinfoptr_1=rpos1->second;
			for(rpos2=rpos1;rpos2!=reslist->resmap.end();rpos2++){
				resinfoptr_2=rpos2->second;
					//create a new CInelasticInfo object from the two current resonance info
				temp = new CInelasticInfo(resinfoptr_1, resinfoptr_2, 0);
				SortedAdd(ThermalArray[abs(temp->net_b)][abs(temp->net_q)][abs(temp->net_s)][Misc::Sign(temp->net_b)][Misc::Sign(temp->net_q)][Misc::Sign(temp->net_s)], *temp);
				delete temp;
			}
		}
		if(UseInelasticArray){
			for(rpos1=reslist->resmap.begin();rpos1!=reslist->resmap.end();rpos1++){
				resinfoptr_1=rpos1->second;
				for(rpos2=reslist->resmap.begin();rpos2!=reslist->resmap.end();rpos2++){
					resinfoptr_2=rpos2->second;
					temp = new CInelasticInfo(resinfoptr_1, resinfoptr_2, 0);
					Th_iter = ThermalArray[abs(temp->net_b)][abs(temp->net_q)][abs(temp->net_s)][Misc::Sign(temp->net_b)][Misc::Sign(temp->net_q)][Misc::Sign(temp->net_s)].begin();
					while(Th_iter != ThermalArray[abs(temp->net_b)][abs(temp->net_q)][abs(temp->net_s)][Misc::Sign(temp->net_b)][Misc::Sign(temp->net_q)][Misc::Sign(temp->net_s)].end()){
						ires1 = Th_iter->resinfo_1->ires;
						ires2 = Th_iter->resinfo_2->ires;
						if(AddToArrayCheck(*(Th_iter->resinfo_1), *(Th_iter->resinfo_2), *resinfoptr_1, *resinfoptr_2)){
							SortedAdd(InelasticArray[min(ires1, ires2)][max(ires1, ires2)], *temp);
							foo++;
							if(foo >= 5000000000){
								CLog::Fatal("CInelastic, Error, too many resonances are being read in\n");
							}
						}
						Th_iter++;

					}
					delete temp;
				}
			}
		}
		if(UseFile){
			inelasticfile.open(filename.c_str(), fstream::out);
			if(inelasticfile){
				//inelasticfile << boltzmann->run_name << endl;
				inelasticfile << "#THERMAL" << endl;
				for(int i = 0; i<3; i++){
					for(int j = 0; j<5; j++){
						for(int k = 0; k<7; k++){
							for(int pm1 = 0; pm1 <2; pm1++){
								for(int pm2 = 0; pm2<2; pm2++){
									for(int pm3 = 0; pm3<2; pm3++){
										inelasticfile << i <<"\t"<< j << "\t"<< k << "\t"<< pm1 << "\t" << pm2 <<"\t"<< pm3 << "\t" << ThermalArray[i][j][k][pm1][pm2][pm3].size() << endl;
										Th_iter = ThermalArray[i][j][k][pm1][pm2][pm3].begin();
										while(Th_iter != ThermalArray[i][j][k][pm1][pm2][pm3].end()){
											inelasticfile << "\t" << Th_iter->resinfo_1->pid << "\t" << Th_iter->resinfo_2->pid << endl;
											Th_iter++;
										}
										inelasticfile << endl;
									}
								}
							}
						}
					}
				}
				if(UseInelasticArray){
					inelasticfile << "#INELASTIC_ARRAY" << endl;
					for(int i = 1; i <= nres; i++){
						for(int j = 1; j <= nres; j++){
							list<CInelasticInfo> templist = InelasticArray[i][j];
							if(templist.size() > 0){
								inelasticfile << i << "\t" << j << "\t" << templist.size() << endl;
								Th_iter = InelasticArray[i][j].begin();
								while(Th_iter != InelasticArray[i][j].end()){
									inelasticfile << "\t" << Th_iter->resinfo_1->ires << "\t" << Th_iter->resinfo_2->ires << endl;
								}
							}
						}
					}
				}
			}else{
				CLog::Fatal("Unable to open inelastic file "+filename+"\n");
			}
			inelasticfile.close();	
		}
	}
}

// adds CInelasticInfo inelastic_in to list list_in, assuming that list_in is sorted by energy
void CInelasticList::SortedAdd(list<CInelasticInfo> &list_in, CInelasticInfo inelastic_in){
	list<CInelasticInfo>::iterator CLiterator = list_in.begin();

	while(CLiterator->min_mass < inelastic_in.min_mass && CLiterator != list_in.end() ){
		CLiterator++;
	}
	list_in.insert(CLiterator, inelastic_in);
}

/*
Method to check whether or not a given inelastic scattering is valid. Allows for more complex checks.
Assumes that res1 and res2 are the incoming particles, while res3 and res4 are outgoing particles.
*/
bool CInelasticList::AddToArrayCheck(CresInfo res1, CresInfo res2, CresInfo res3, CresInfo res4){
	bool add = true;
	//check for charge, baryon number, and strangness conservation
	if(res1.charge + res2.charge - res3.charge - res4.charge != 0){
		add = false;
	}
	if(res1.baryon + res2.baryon - res3.baryon - res4.baryon != 0){
		add = false;
	}
	if(res1.strange + res2.strange - res3.strange - res4.strange != 0){
		add = false;
	}

	//check that the incoming and outgoing resonances aren't the same
	if((res1.ires == res3.ires && res2.ires == res4.ires)||(res1.ires == res4.ires && res2.ires == res3.ires)){
		add = false;
	}

	return add;
}
