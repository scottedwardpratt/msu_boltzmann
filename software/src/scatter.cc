#include "msu_commonutils/misc.h"
#include "msu_commonutils/randy.h"
#include "msu_commonutils/constants.h"
#include "msu_boltzmann/msu_boltzmann.h"
#include "msu_boltzmann/msupart.h"
#include "msu_sampler/resonances.h"

//using namespace Misc;

/*
InelasticScatter: Method to implement inelastic scattering. Right now, it magically
converts the two incoming particles to outgoing particles and then scatters them
elastically.
*/

void CMSU_Boltzmann::InelasticScatter(CMSUPart *part1,CMSUPart *part2,CMSUPart *part3,CMSUPart *part4,CInelasticInfo inelinfo){
	double m0squared,m3squared,m4squared,m0,mtot,cthet,sthet,phi,q;
	FourVector *p3=&part3->p,*p4=&part4->p;
	FourVector p0tot,u,pprime;
	const FourVector g={1,-1,-1,-1};
	int alpha;

	part3->resinfo = inelinfo.resinfo_1;
	part4->resinfo = inelinfo.resinfo_2;

	m0squared=0;
	for(int b = 0; b< 4; b++){
		p0tot[b] = part1->p[b]+part2->p[b];
		m0squared+=g[b]*p0tot[b]*p0tot[b];
	}
	m3squared=pow(part3->resinfo->mass,2);
	m4squared=pow(part4->resinfo->mass,2);
	mtot=sqrt(m3squared)+sqrt(m4squared);

	if(mtot*mtot<m0squared){
		//pick the angles
		cthet=1.0-2.0*randy->ran();
		sthet=sqrt(1.0-cthet*cthet);
		phi=2.0*PI*randy->ran();
		q=sqrt(Misc::triangle2(m0squared,m3squared,m4squared));
		(*p3)[3]=q*cthet;
		(*p3)[1]=q*sthet*cos(phi);
		(*p3)[2]=q*sthet*sin(phi);
		(*p3)[0]=sqrt(m3squared+(*p3)[1]*(*p3)[1]+(*p3)[2]*(*p3)[2]+(*p3)[3]*(*p3)[3]);
		(*p4)[3]=-(*p3)[3];
		(*p4)[2]=-(*p3)[2];
		(*p4)[1]=-(*p3)[1];
		(*p4)[0]=sqrt(m4squared+(*p4)[1]*(*p4)[1]+(*p4)[2]*(*p4)[2]+(*p4)[3]*(*p4)[3]);
		/* Boost the new particles */
		m0=sqrt(m0squared);
		for(alpha=0;alpha<4;alpha++)
			u[alpha]=p0tot[alpha]/m0;
		Misc::lorentz(u,*p3,pprime);
		for(alpha=0;alpha<4;alpha++)
			(*p3)[alpha]=pprime[alpha];
		part3->msquared=m3squared;
		part3->SetY();
		Misc::lorentz(u,*p4,pprime);
		for(alpha=0;alpha<4;alpha++)
			(*p4)[alpha]=pprime[alpha];
		part4->msquared=m4squared;
		part4->SetY();
	}
}

/*
Scatter: Elastically scatters two CMSUPart objects with s-wave angular distribution
*/

void CMSU_Boltzmann::Scatter(CMSUPart *part1,CMSUPart *part2,CMSUPart *part3,CMSUPart *part4){
	double roots=0.0;
	double ctheta,stheta,phi,qmag;
	FourVector ptot,u,q,qprime;
	const FourVector g={1,-1,-1,-1};
	FourVector *p1=&part1->p,*p2=&part2->p,*p3=&part3->p,*p4=&part4->p;
	int alpha;
	if(part3!=part1)
		part3->resinfo=part1->resinfo;
	if(part4!=part2)
		part4->resinfo=part2->resinfo;

	for(alpha=0;alpha<4;alpha++){
		ptot[alpha]=(*p1)[alpha]+(*p2)[alpha];
		q[alpha]=0.5*((*p1)[alpha]-(*p2)[alpha]);
		roots+=g[alpha]*ptot[alpha]*ptot[alpha];
	}
	roots=sqrt(roots);
	for(alpha=0;alpha<4;alpha++)
		u[alpha]=ptot[alpha]/roots;
	Misc::BoostToCM(u,q,qprime);
	qmag=sqrt(qprime[1]*qprime[1]+qprime[2]*qprime[2]+qprime[3]*qprime[3]);

	ctheta=1.0-2.0*randy->ran();
	phi=2.0*PI*randy->ran();
	qprime[3]=qmag*ctheta;
	stheta=sqrt(1.0-ctheta*ctheta);
	qprime[1]=qmag*stheta*cos(phi);
	qprime[2]=qmag*stheta*sin(phi);

	qprime[0]=sqrt(part1->msquared+qmag*qmag);
	Misc::Boost(u,qprime,*p3);
	for(alpha=1;alpha<4;alpha++)
		qprime[alpha]=-qprime[alpha];
	qprime[0]=sqrt(part2->msquared+qmag*qmag);
	Misc::Boost(u,qprime,*p4);
	part3->msquared=part1->msquared;
	part4->msquared=part2->msquared;
	part3->Setp0();
	part4->Setp0();
	part3->SetY();
	part4->SetY();
}

/*
If resonance can decay, A->B+C, then if B+C collide, this will merge B+C->A
*/
bool CMSU_Boltzmann::Merge(CMSUPart *part1,CMSUPart *part2,CMSUPart *part3,CresInfo *resinfo){
	int alpha;
	double ptot[4],s;
	FourVector *p1=&part1->p,*p2=&part2->p;
	bool success;

	for(alpha=0;alpha<4;alpha++)
		ptot[alpha]=(*p1)[alpha]+(*p2)[alpha];
	s=ptot[0]*ptot[0]-ptot[1]*ptot[1]-ptot[2]*ptot[2]-ptot[3]*ptot[3];
	if(s>resinfo->minmass*resinfo->minmass){
		for(alpha=0;alpha<4;alpha++)
			part3->p[alpha]=ptot[alpha];
		part3->msquared=s;
		part3->resinfo=resinfo;
		part3->SetY();
		success=true;
		part3->balanceID=-1;
		//part3->eta0=part1->eta0;
		part3->phi0=part1->phi0;
	}
	else{
		success=false;
	}
	return success;
}
