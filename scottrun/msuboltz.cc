#include "msu_sampler/master.h"
#include "msu_commonutils/qualifier.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include <cstring>
using namespace std;

int main(){
	CparameterMap parmap;
	parmap.ReadParsFromFile("parameters_sampler.txt");
	string run_name="default_0";
	printf("check 000000\n");

	CresList *reslist=new CresList(&parmap);
	printf("check 00000\n");
	CmasterSampler::meanfield=new CmeanField_Simple(&parmap);
		printf("check 0000\n");
	CmasterSampler ms(&parmap);
	CpartList *pl=new CpartList(&parmap,reslist);
		printf("check 000\n");

	ms.partlist=pl;
	printf("check 00\n");
	ms.randy->reset(time(NULL));
	printf("check 0\n");
	ms.ReadHyper_OSU_2D();

	printf("check a\n");
	CMSU_Boltzmann *msuboltz=new CMSU_Boltzmann(run_name,reslist);
	printf("check b\n");
	msuboltz->InitCascade();
	printf("check c\n");
	
	long long int nparts=0;
	int	nevents=parmap.getI("SAMPLER_NEVENTS",10);
	nevents=1;

	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.txt");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		for(int ievent=0;ievent<nevents;ievent++){
			printf("check d\n");
			nparts+=ms.MakeEvent();
			printf("check e\n");
			msuboltz->InputPartList(pl);
			pl->Clear();
			printf("ievent=%lld nparts/event=%g\n",ms.NEVENTS,double(nparts)/double(ms.NEVENTS));
		}
	}

	printf("YIPPEE!!!!! We made it all the way through!\n");
	return 0;
}
