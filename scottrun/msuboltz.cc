#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include <cstring>
using namespace std;

int main(){
	CparameterMap parmap;
	string run_name="default_0";

	string filename="model_output/fixed_parameters.txt";
	parmap.ReadParsFromFile(filename);
	filename="model_output/"+run_name+"/parameters.txt";
	parmap.ReadParsFromFile(filename);

	CmasterSampler ms(&parmap);
	CpartList *pl=new CpartList(&parmap,ms.reslist);

	ms.partlist=pl;
	ms.randy->reset(time(NULL));
	ms.ReadHyper_OSU_2D();

	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_name,&parmap,ms.reslist);
	msuboltz->InitCascade();
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	nevents=parmap.getI("MSU_BOLTZMANN_NEVENTSMAX",10);

	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		for(int ievent=0;ievent<nevents;ievent++){
			msuboltz->Reset();
			nparts+=ms.MakeEvent();
			msuboltz->InputPartList(pl);
			pl->Clear();
			msuboltz->PerformAllActions();
			printf("ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
		}
	}

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
