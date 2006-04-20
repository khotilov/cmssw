#include "SimG4Core/CustomPhysics/interface/CustomPDGParser.h"

bool CustomPDGParser::s_isRHadron(int pdg) 
{
    int pdgAbs=abs(pdg);
    return ((pdgAbs % 100000 / 10000 == 9) ||  
	    (pdgAbs % 10000 / 1000 == 9) || s_isRGlueball(pdg));
}

bool CustomPDGParser::s_isSLepton(int pdg)
{
    int pdgAbs=abs(pdg);
    return (pdgAbs / 100 % 10000 == 0 && pdgAbs / 10 % 10 == 1);
}

bool CustomPDGParser::s_isRBaryon(int pdg)
{
    int pdgAbs=abs(pdg);
    return (pdgAbs % 100000 / 10000 == 9);
}

bool CustomPDGParser::s_isRGlueball(int pdg)
{
    int pdgAbs=abs(pdg);
    return (pdgAbs % 100000 / 10 == 99);

}

bool CustomPDGParser::s_isRMeson(int pdg)
{
    int pdgAbs=abs(pdg);
    return (pdgAbs % 10000 / 1000 == 9);

}

double CustomPDGParser::s_charge(int pdg)
{
    float charge=0,sign=1;
    int pdgAbs=abs(pdg);
    if(pdg < 0 ) sign=-1;

    if (s_isSLepton(pdg))     //Sleptons
    {
	if( pdgAbs %2 == 0) return 0;
	else return -sign;
    }

    if (s_isRMeson(pdg))
    {
        std::vector<int> quarks = s_containedQuarks(pdg);
        if (quarks[1] % 2 == 0 && quarks[0] % 2 == 1) charge=1;
        if (quarks[1] % 2 == 1 && quarks[0] % 2 == 0 ) charge=-1;
        charge*=sign;       
	return charge;
    }

    if (s_isRBaryon(pdg))
    {
	int baryon = s_containedQuarksCode(pdg);
	for (int q=1; q< 1000; q*=10)
	{
	    if(baryon / q % 2 == 0) charge+=2; else charge -=1; 
	}
        charge/=3;
	charge*=sign;
	return charge;
    }

    return 0; 
}

double CustomPDGParser::s_spin(int pdg)
{
    int pdgAbs=abs(pdg);
    return pdgAbs % 10;    
}

std::vector<int> CustomPDGParser::s_containedQuarks(int pdg)
{
    std::vector<int> quarks;
    for (int i=s_containedQuarksCode(pdg); i > 0; i /= 10)
    {
	quarks.push_back(i % 10);
    }
    return quarks; 
}

int CustomPDGParser::s_containedQuarksCode(int pdg)
{
    int pdgAbs=abs(pdg);
    if (s_isRBaryon(pdg))
	return pdgAbs / 10 % 1000;

    if(s_isRMeson(pdg))
	return pdgAbs / 10 % 100;

    return 0;
}
