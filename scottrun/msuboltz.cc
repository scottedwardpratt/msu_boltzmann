#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include <cstring>
using namespace std;

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters_sampler.txt");
	string run_name="default_0";

	CresList *reslist=new CresList(&parmap);
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
	CmasterSampler ms(&parmap);
	CpartList *pl=new CpartList(&parmap,reslist);

	ms.partlist=pl;
	ms.randy->reset(time(NULL));
	ms.ReadHyper_OSU_2D();

	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_name);
	msuboltz->InitCascade();
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	nevents=1;

	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		for(int ievent=0;ievent<nevents;ievent++){
			nparts+=ms.MakeEvent();
			msuboltz->InputPartList(pl);
			pl->Clear();
			printf("ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
		}
	}

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
