#ifndef PARTICLEDEPOSIT_HH_
#define PARTICLEDEPOSIT_HH_

#include <vector>

#ifdef __MAKECINT__
#pragma link C++ class vector<minmiser::ParticleDeposit>
#endif

#include "RecoParticleFlow/PFClusterTools/interface/Deposition.h"
#include "RecoParticleFlow/PFClusterTools/interface/DetectorElement.h"


#include <iostream>

/**
 \class ParticleDeposit 
 \brief An object that encapsualtes energy depositions (real, MC or otherwise) made by  particle in N DetectorElement objects.

 \author Jamie Ballin
 \date   April 2008
 */
namespace pftools {

class ParticleDeposit {
public:

	ParticleDeposit(double truthEnergy = -1.0, double eta = 0, double phi = 0);
	virtual ~ParticleDeposit();

	void addRecDeposition(Deposition rec);
	void addTruthDeposition(Deposition truth);

	const std::vector<Deposition>& getRecDepositions() const;
	std::vector<Deposition> getTruthDepositions() const;

	/*
	 * Returns the overall MC particle energy.
	 */
	inline double getTruthEnergy() const {
		return myTruthEnergy;
	}

	/* 
	 * Returns the detected energy from this detector element, including calibration.
	 */
	double getRecEnergy(const DetectorElementPtr de) const;

	double getRecEnergy() const;

	double getEnergyResolution() const;

	/*
	 * Returns the raw MC energy input into this detector element.
	 */
	double getTruthEnergy(const DetectorElementPtr de) const;

	inline unsigned getId() const {
		return myId;
	}

	inline double getEta() const {
		return myEta;
	}

	inline double getPhi() const {
		return myPhi;
	}

	friend std::ostream& operator<<(std::ostream& s, const ParticleDeposit& p);

private:
	static unsigned count;
	//ParticleDeposit(const ParticleDeposit& pd);
	std::vector<Deposition> myRecDepositions;
	std::vector<Deposition> myTruthDepositions;
	unsigned myId;

	double myTruthEnergy;
	double myEta;
	double myPhi;

	/*
	 * For general ROOT dictionary building happiness!
	 */
	std::vector<ParticleDeposit*> pdps_;
	//	std::vector<ParticleDeposit> pds_;
};


typedef boost::shared_ptr<ParticleDeposit> ParticleDepositPtr;

}

#endif /*PARTICLEDEPOSIT_HH_*/
