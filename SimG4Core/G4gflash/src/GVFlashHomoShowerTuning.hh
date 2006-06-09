//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef GVFlashHomoShowerTuning_hh
#define GVFlashHomoShowerTuning_hh

// J.P. Wellisch, Oct. 2004

class GVFlashHomoShowerTuning
{
	public:
	// Definitions:
	// <t>: shower center of gravity
	// T: Depth at shower maximum
	// Ec: Critical energy
	// y = E/Ec
	// X0: Radiation length
	
	// Homogeneous media:
	// Avarage shower profile
	// (1/E)(dE(t)/dt) = f(t) = (beta*t)**(alpha-1)*beta*std::exp(-beta*t)/Gamma(alpha)
	// where Gamma: The Gamma function
	//
	// <t> = alpha/beta
	// T = (alpha-1)/beta
	// and
	// T = ln(y) + t1
	// alpha = a1+(a2+a3/Z)ln(y)
	virtual G4double ParAveT1(){ return -0.858;} // t1
	virtual G4double ParAveA1(){ return 0.21;  } // a1
	virtual G4double ParAveA2(){ return 0.492; } // a2
	virtual G4double ParAveA3(){ return 2.38;  } // a3
	
	// std::sqrt(var(ln(T))) = 1/(t+t2*ln(y))
	virtual G4double ParSigLogT1(){ return -1.4;} // t1
	virtual G4double ParSigLogT2(){ return 1.26;} // t2
	
	// std::sqrt(var(ln(alpha))) = 1/(a1+a2*ln(y))
	virtual G4double ParSigLogA1(){ return -0.58;} // a1
	virtual G4double ParSigLogA2(){ return 0.86; } // a2
	
	// Correlation(ln(T),ln(alpha))=r1+r2*ln(y)
	virtual G4double ParRho1(){ return 0.705; } // r1
	virtual G4double ParRho2(){ return -0.023;} // r2
	
	// Radial profiles
	// f(r) := (1/dE(t))(dE(t,r)/dr)
	// Ansatz:
	// f(r) = p(2*r*Rc**2)/(r**2+Rc**2)**2+(1-p)*(2*r*Rt**2)/(r**2+Rt**2)**2, 0<p<1
	// 
	// Rc (t/T)= z1 +z2*t/T
	// z1 = c1+c2*ln(E/GeV)
	// z2 = c3+c4*Z
	virtual G4double ParRC1(){ return 0.0251;   } // c1
	virtual G4double ParRC2(){ return 0.00319;  } // c2
	virtual G4double ParRC3(){ return 0.1162;   } // c3
	virtual G4double ParRC4(){ return -0.000381;} // c4
	
	// Rt (t/T)= k1*(std::exp(k3*(t/T-k2))+std::exp(k4*(t/T-k2)))
	// k1 = t1+t2*Z
	// k4 = t5+t6*ln(E/GeV)
	virtual G4double ParRT1(){ return 0.659;   } // t1
	virtual G4double ParRT2(){ return -0.00309;} // t2
	virtual G4double ParRT3(){ return 0.645;   } // k2
	virtual G4double ParRT4(){ return -2.59;   } // k3
	virtual G4double ParRT5(){ return 0.3585;  } // t5
	virtual G4double ParRT6(){ return 0.0412;  } // t6
	
	// p(t/T) = p1*std::exp((p2-t/T)/p3 - std::exp((p2-t/T)/p3))
	// p1 = c1+c2*Z
	// p2 = c3+c4*Z
	// p3 = c5 + c6*ln(E/GeV)
	virtual G4double ParWC1(){ return 2.632;   } // c1
	virtual G4double ParWC2(){ return -0.00094;} // c2
	virtual G4double ParWC3(){ return 0.401;   } // c3
	virtual G4double ParWC4(){ return 0.00187; } // c4
	virtual G4double ParWC5(){ return 1.313;   } // c5
	virtual G4double ParWC6(){ return -0.0686; } // c6
	
	// Fluctuations on radial profiles through number of spots
	// The total number of spots needed for a shower is
	// Ns = n1*ln(Z)(E/GeV)**n2
	virtual G4double ParSpotN1(){ return 93.;  } // n1
	virtual G4double ParSpotN2(){ return 0.876;} // n2
	
	// The number of spots per longitudinal interval is:
	// (1/Ns)(dNs(t)/dt) = f(t) = (beta*t)**(alpha-1)*beta*std::exp(-beta*t)/Gamma(alpha)
	// <t> = alpha_s/beta_s
	// Ts = (alpha_s-1)/beta_s
	// and
	// Ts = T*(t1+t2*Z)
	// alpha_s = alpha*(a1+a2*Z)
	virtual G4double ParSpotT1(){ return 0.698;  } // t1
	virtual G4double ParSpotT2(){ return 0.00212;} // t2
	
	virtual G4double ParSpotA1(){ return 0.639;  } //a1
	virtual G4double ParSpotA2(){ return 0.00334;} //a2
	
};

#endif


