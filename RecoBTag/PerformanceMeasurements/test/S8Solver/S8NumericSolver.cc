/////////////////////////////////////////////////////////////////////////
//
//  System is a method to determine selection effiencies and
//  signal/background fractions on data with minimum monte-carlo input
//
//  It requires :
//     - 1 data sample containing 2 signal components
//     - 3 selection criteria (cuts on sensitive variable) with
//       various efficiencies on the 2 signals.
//
//  The solved system is :
//   / 1    = na + nb
//  |  q1   = ea1*na + eb1*nb
//  |  q2   = ea2*na + eb2*nb
//  /  q3   = ea3*na + eb3*nb
//  \\  q12  = ca12*ea1*ea2*na + cb12*eb1*eb2*nb
//  |  q23  = ca23*ea2*ea3*na + cb23*eb2*eb3*nb
//  |  q31  = ca31*ea3*ea1*na + cb31*eb3*eb1*nb
//   \\ q123 = ca12*ca23*c31*ea1*ea2*ea3*na + cb12*cb23*cb31**eb2*eb3*nb
//
//   Input :
//     qxxx        : fractions of data surviving selection(s) xxx
//     caxy (cbxy) : correlation factor between cuts x and y on signal a (b)
//
//   Output :
//     na (nb)     : fraction signal a (b)
//     eax (ebx)   : efficiency of cut x on signal a (b)
//
/////////////////////////////////////////////////////////////////////////
#include "S8NumericSolver.h"


S8NumericSolver::S8NumericSolver() : TNamed()
{
// Generic constructor
this->Reset();
}


S8NumericSolver::S8NumericSolver(std::string name) : TNamed(name,name)
{
this->Reset();
}

S8NumericSolver::S8NumericSolver(std::string name, Double_t n, Double_t n1, Double_t n2 , Double_t n3,
                             Double_t n12, Double_t n23, Double_t n31, Double_t n123) : TNamed(name,name)
{
// Constructor with input initialization :
// n    : total numer of events in data sample
// nx   : number of events passing cut x
// nxy  : number of events passing cuts x and y
// n123 : number of events passing all 3 cuts.
  this->Reset();
  this->SetInput(n,n1,n2,n3,n12,n23,n31,n123);
}

S8NumericSolver::S8NumericSolver(std::string name, Double_t* n) : TNamed(name,name)
{
this->Reset();
this->SetInput(n[0],n[1],n[2],n[3],n[4],n[5],n[6],n[7]);
}


void S8NumericSolver::Reset()
{
// Reseting inputs, correlations, results and errors
kError = 0;
fRndm = new TRandom3(12345);

for(Int_t i=0 ; i<8 ; ++i)
  {
  fInput[i]  = 0.  ;
  fIndep[i]  = 0.  ;

  fErrorinf_Stat[i] = 0. ;
  fErrorsup_Stat[i] = 0.;
  fErrorinf_Syst[i] = 0. ;
  fErrorsup_Syst[i] = 0. ;

  if( i<6 )
    {
    fCorr[i/3][i%3]    = 1. ;
    fCorrerr[i/3][i%3] = 1. ;
    kc[i/3][i%3]= 1.;
    if(i<3)
      {
      q[i]=-1.;
      Q[i]=-1.;
      }
    }

  }
QQQ= -1. ;
fAsym[0]= 0;
fAsym[1]= 1;
 fNpt = 20000;//5000;
fIter = 10000;
 fAveRes = 0.;
 fAveResSetup = false;
 fverbose = false;
 
 fForceSol = false;
 fpickSol = 0;

}

void S8NumericSolver::SetInput(Double_t n, Double_t n1, Double_t n2 , Double_t n3,
                             Double_t n12, Double_t n23, Double_t n31, Double_t n123)
{
// Set the 8 input values - cf constructor
  fInput[0] = n  ;
  fInput[1] = n1  ;
  fInput[2] = n2  ;
  fInput[3] = n3  ;
  fInput[4] = n12  ;
  fInput[5] = n23  ;
  fInput[6] = n31  ;
  fInput[7] = n123  ;

  fIndep[0] = sqrt(n123)  ;
  fIndep[1] = sqrt(n23 - fIndep[0])  ;
  fIndep[2] = sqrt(n31 - fIndep[0])  ;
  fIndep[3] = sqrt(n3 - fIndep[0] - fIndep[1] - fIndep[2]) ;
  fIndep[4] = sqrt(n12 - fIndep[0])  ;
  fIndep[5] = sqrt(n2 - n12 - fIndep[1])  ;
  fIndep[6] = sqrt(n1 - n12 - fIndep[2])  ;
  fIndep[7] = sqrt(n - fIndep[6] - fIndep[5] - fIndep[4] - n3) ;
}

void S8NumericSolver::SetCorr(Double_t c12a, Double_t c23a, Double_t c31a,
                            Double_t c12b, Double_t c23b, Double_t c31b)
{
// Set the correlation factors
// 3 first parameters : correlation for signal a.
// 3 last parameters  : correlation for signal b.

  fCorr[0][0] = c12a  ;
  fCorr[1][0] = c12b  ;
  fCorr[0][1] = c23a  ;
  fCorr[1][1] = c23b  ;
  fCorr[0][2] = c31a  ;
  fCorr[1][2] = c31b  ;
}

void S8NumericSolver::SetCorrError(Double_t c12a, Double_t c23a, Double_t c31a,
                            Double_t c12b, Double_t c23b, Double_t c31b)
{
// Set the errors on correlation factors
// Same order than SetCorr
  fCorrerr[0][0] = c12a  ;
  fCorrerr[1][0] = c12b  ;
  fCorrerr[0][1] = c23a  ;
  fCorrerr[1][1] = c23b  ;
  fCorrerr[0][2] = c31a  ;
  fCorrerr[1][2] = c31b  ;
}


void S8NumericSolver::MakeSystem(Double_t *shift)
{
Double_t N = fInput[0]+shift[0];
for(Int_t i=0;i<3;++i)
  {
  q[i] = (fInput[i+1]+shift[i+1])/N ;
  Q[i] = (fInput[i+4]+shift[i+4])/N ;
  }
QQQ = fInput[7]/N;
for(Int_t i=0;i<6;++i)
  {
  kc[i/3][i%3] = fCorr[i/3][i%3] + shift[i+8];
  }
}

Int_t S8NumericSolver::Solve()
{
Double_t shift[14]={0};
Double_t res[5]={0};
MakeSystem(shift);

 fverbose = true; 
if(!FindSolution(res,SolveSystem(res))) {
	std::cout << "[S8Numeric SOLUTIONS NOT FOUND, leaving ...." << std::endl;
	return 0;
}
 std::cout << "[S8Numeric] got solutions, now compute errors" << std::endl;
 fverbose = false;
 ComputeErrors();
return 1;
}

Int_t S8NumericSolver::FindSolution(Double_t* res, int n)
{

	int npositiveSols = 0;
	int nphysicalSols = 0;
	int thesols = -1;
	int thesols2 = -1;
	double deltares = 9999999;

	if (fForceSol) {
	  if (fverbose) std::cout << "[S8Numeric] Solution chosen manually, Force to be solution # " << fpickSol<< std::endl;
	  fNb = res[fpickSol -1 ];
	  for(int i=0; i<8;++i) fMapResult[i] =  fResult[i] = E(i%2,i/2);
	  return 1;
	}

	if (fverbose) std::cout << "[S8Numeric] now print all solutions:"<< std::endl;
	for(int j=1; j<= n; ++j) {
		fNb = res[j-1];
		if (fverbose) std::cout << "solution n= " << j << std::endl;
		double tmpsol[8];
		double totprod=1;
		for(int i=0; i<8;++i) {
			tmpsol[i] = E(i%2,i/2);
			if (fverbose) std::cout << " result i=" << i << " " <<  tmpsol[i] << std::endl;
			totprod=totprod*tmpsol[i];
		}
		if (totprod>=0) {
			npositiveSols++;
			int tmpcounter = 0;
			for(int ii=2; ii<=5;++ii) {
				
				if (tmpsol[ii]>=0 && tmpsol[ii]<=1.) tmpcounter++;
			}
			if (tmpcounter==4) nphysicalSols++;
			double tmpdeltares = fabs(fAveRes - res[j-1]);
			if ( nphysicalSols>0 ) {
				if (fAveResSetup && (deltares>tmpdeltares) ) {thesols = j-1; deltares = tmpdeltares;}
				else if (!fAveResSetup && tmpcounter==4 ) { 
				thesols = j-1;}
				else if (!fAveResSetup && nphysicalSols>1 && tmpcounter==4 ) {
					if ( tmpsol[4]>0 && tmpsol[4]>tmpsol[5] && tmpsol[2]>tmpsol[3] )				thesols = j-1;
				}
			}
		}
		if ( tmpsol[4]>0 && tmpsol[4]>tmpsol[5] && tmpsol[2]>tmpsol[3] ) thesols2= j-1;
		
		if (fverbose) {
			std::cout << "\n";
			std::cout << "# physical solutions = " << nphysicalSols << std::endl;
		}
	}
		
	if (nphysicalSols==0) {
		if ( thesols2 != -1 ) {
			fNb = res[thesols2];
			for(int i=0; i<8;++i) fMapResult[i] =  fResult[i] = E(i%2,i/2);
		}
		else return 0;
	}
	else {
		std::cout << "the sols = " << thesols << std::endl;
		fNb = res[thesols];
		for(int i=0; i<8;++i) fMapResult[i] =  fResult[i] = E(i%2,i/2);
		return 1;
	}
/*
	
if(n<=0)
  {
//  printf("\nNo physics solution.\n");
  return 0;
  }
else if(n==1) fNb = res[0];
else if(n==2)
  {
  Double_t a=0,b=0;
  Bool_t bb[2];
  for(Int_t i=0;i<2;++i)
    {
    fNb = res[i];
    Int_t u = (2*i-1)*fAsym[1];
    a = E(0,fAsym[0])*u;
    b = E(1,fAsym[0])*u;
    bb[i] = (a<b);
    }
  if     ( bb[0] &&  bb[1]) fNb = res[0];
  else if(!bb[0] && !bb[1]) fNb = res[1];
  else
    {
    printf("\nCannot solve ambiguity... Try SetInitialOrder on another parameter\n");
    return 0;
    }
  }
else if(n>2)
  {
  printf("\nToo many solutions. What's going on here ?\n");
  return 0;
  }

*/
	
//for(int i=0; i<8;++i) fResult[i] = E(i%2,i/2);
//return 1;
}


void S8NumericSolver::ComputeErrors()
{
printf("Start computing errors -> mode : %d\n",kError);
if(kError==0) return;

Double_t central[8]={0};
Double_t res[5];
for(int i=0;i<8;++i)
  {
	  central[i] = fMapResult[i];//fResult[i];
  fErrorsup_Stat[i] = 0.;
  fErrorinf_Stat[i] = 0.;
  fErrorsup_Syst[i] = 0.;
  fErrorinf_Syst[i] = 0.;
  fMapErrorInf_Stat[i] = 0;
  fMapErrorSup_Stat[i] = 0;
  
  }
Int_t Nsup = 0;
Int_t Ninf = 0;
Double_t shift[14] = {0.};
if(kError==1 || kError==2) // Stat error
  {
  printf("   Computing stat errors\n");
  Ninf = Nsup = 0;
  for(int i=8;i<14;++i) shift[i] = 0.;
  for(int i=0;i<fIter;++i)
    {
    Double_t w[8] = {0};

    for(int j=0;j<8;++j) w[j] = fRndm->Gaus(0.,1.)*fIndep[j];
    shift[0] = w[0]+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7];
    shift[1] = w[0] + w[2] + w[4] + w[6];
    shift[2] = w[0] + w[1] + w[4] + w[5];
    shift[3] = w[0] + w[1] + w[2] + w[3];
    shift[4] = w[0] + w[4];
    shift[5] = w[0] + w[1];
    shift[6] = w[0] + w[2];
    shift[7] = w[0] ;
    MakeSystem(shift);

    if(!FindSolution(res,SolveSystem(res))) continue;
    for(int j=0;j<8;++j)
      {
		  //if(fResult[j]<central[j] && fResult[j]>0) {fErrorinf_Stat[j]+= pow(fResult[j]-central[j],2); Ninf++;}
		  //if(fResult[j]>central[j] && fResult[j]<1) {fErrorsup_Stat[j]+= pow(fResult[j]-central[j],2); Nsup++;}
	  if(fMapResult[j]<central[j] && fMapResult[j]>0) {fMapErrorInf_Stat[j]+= pow(fMapResult[j]-central[j],2); Ninf++;}
      if(fMapResult[j]>central[j] && fMapResult[j]<1) {fMapErrorSup_Stat[j]+= pow(fMapResult[j]-central[j],2); Nsup++;}
      }
    }
  for(int i=0;i<8;++i)
    {
		//fErrorinf_Stat[i]/=Ninf;
		//fErrorsup_Stat[i]/=Nsup;
	fMapErrorInf_Stat[i]/=Ninf;
    fMapErrorSup_Stat[i]/=Nsup;  
    }
  }
if(kError==1 || kError==3) // Syst error
  {
  printf("   Computing syst errors\n");
  Ninf = Nsup = 0;
  for(int i=0;i<8;++i) shift[i] = 0.;
  for(int i=0;i<fIter;++i)
    {
    for(Int_t j=0;j<6;++j) shift[j+8] = fRndm->Gaus(0.,1.)*fCorrerr[j/3][j%3];
    MakeSystem(shift);
    if(!FindSolution(res,SolveSystem(res))) continue;
    for(int j=0;j<8;++j)
      {
		  //if(fResult[j]<central[j] && fResult[j]>0) {fErrorinf_Syst[j]+= pow(fResult[j]-central[j],2); Ninf++;}
		  //if(fResult[j]>central[j] && fResult[j]<1) {fErrorsup_Syst[j]+= pow(fResult[j]-central[j],2); Nsup++;}
	  if(fMapResult[j]<central[j] && fMapResult[j]>0) {fErrorinf_Stat[j]+= pow(fMapResult[j]-central[j],2); Ninf++;}
      if(fMapResult[j]>central[j] && fMapResult[j]<1) {fErrorsup_Stat[j]+= pow(fMapResult[j]-central[j],2); Nsup++;}
      }
    }
  for(int i=0;i<8;++i)
    {
		//fErrorinf_Syst[i]/=Ninf;
		//fErrorsup_Syst[i]/=Nsup;
	fMapErrorInf_Stat[i]/=Ninf;
    fMapErrorSup_Stat[i]/=Nsup;
    }
  }
for(int i=0;i<8;++i) fMapResult[i] = fResult[i] = central[i];
 
}

void S8NumericSolver::SetError( Int_t b )
{
// Error computing mode
// 0 : no errors
// 1 : Stat + Syst
// 2 : Stat only
// 3 : Syst only

kError = b ;
}
Double_t S8NumericSolver::GetValue(Double_t* tab, std::string label)
{
if(label == "na")  return tab[0];
if(label == "nb")  return tab[1];
if(label == "ea1") return tab[2];
if(label == "eb1") return tab[3];
if(label == "ea2") return tab[4];
if(label == "eb2") return tab[5];
if(label == "ea3") return tab[6];
if(label == "eb3") return tab[7];
return 0;
}

Double_t S8NumericSolver::GetResult(Int_t n)
{
// return the system 8 results :
// n = 0 -> na
// n = 1 -> nb
// n = 2 -> ea1
// n = 3 -> eb1
// n = 4 -> ea2
// n = 5 -> eb2
// n = 6 -> ea3
// n = 7 -> eb3

if(n>=0 && n<8) return fResult[n];
return 0;
}

Double_t S8NumericSolver::GetErrorSup(Int_t n ,std::string opt)
{
Double_t err = 0.;
if(n>=0 && n<8 && kError>0)
  {
  if(opt=="All" || opt=="" || opt == "Stat") err += fErrorsup_Stat[n];
  if(opt=="All" || opt=="" || opt == "Syst") err += fErrorsup_Syst[n];
  }
return sqrt(err);
}

Double_t S8NumericSolver::GetErrorInf(Int_t n ,std::string opt)
{
Double_t err = 0.;
if(n>=0 && n<8 && kError>0)
  {
  if(opt=="All" || opt=="" || opt == "Stat") err += fErrorinf_Stat[n];
  if(opt=="All" || opt=="" || opt == "Syst") err += fErrorinf_Syst[n];
  }
return sqrt(err);
}

Double_t S8NumericSolver::GetResult(std::string label)
{
// return the system 8 results :
// label must be one of those
// "na", "nb", "ea1", "eb1", "ea2", "eb2", "ea3", "eb3"
return this->GetValue(fResult, label);
}

Double_t S8NumericSolver::GetErrorSup(std::string label,std::string opt)
{
Double_t err = 0.;
if(opt=="All" || opt=="" || opt == "Stat") err+= this->GetValue(fErrorsup_Stat, label);
if(opt=="All" || opt=="" || opt == "Syst") err+= this->GetValue(fErrorsup_Stat, label);
return sqrt(err);
}

Double_t S8NumericSolver::GetErrorInf(std::string label,std::string opt)
{
Double_t err = 0.;
if(opt=="All" || opt=="" || opt == "Stat") err+= this->GetValue(fErrorinf_Stat, label);
if(opt=="All" || opt=="" || opt == "Syst") err+= this->GetValue(fErrorinf_Stat, label);
return sqrt(err);
}

Double_t S8NumericSolver::GetResultVec(Int_t n)
{
// return the system 8 results :
// n = 0 -> na
// n = 1 -> nb
// n = 2 -> ea1
// n = 3 -> eb1
// n = 4 -> ea2
// n = 5 -> eb2
// n = 6 -> ea3
// n = 7 -> eb3

  //std::cout << "[S8NumericSolver] GetResultVec " << n << " " << fMapResult[n] << std::endl; 
if(n>=0 && n<8) return fMapResult[n];
return 0;
}

Double_t S8NumericSolver::GetErrorSupVec(Int_t n)
{
// return the system 8 results :
// n = 0 -> na
// n = 1 -> nb
// n = 2 -> ea1
// n = 3 -> eb1
// n = 4 -> ea2
// n = 5 -> eb2
// n = 6 -> ea3
// n = 7 -> eb3

  //std::cout << "[S8NumericSolver] GetErrorSupVec " << n << " " << fMapErrorSup_Stat[n] << std::endl; 
if(n>=0 && n<8) return sqrt(fMapErrorSup_Stat[n]);
return 0;
}
Double_t S8NumericSolver::GetErrorInfVec(Int_t n)
{
// return the system 8 results :
// n = 0 -> na
// n = 1 -> nb
// n = 2 -> ea1
// n = 3 -> eb1
// n = 4 -> ea2
// n = 5 -> eb2
// n = 6 -> ea3
// n = 7 -> eb3

  // std::cout << "[S8NumericSolver] GetErrorInfVec " << n << " " << fMapErrorInf_Stat[n] << std::endl; 
if(n>=0 && n<8) return sqrt(fMapErrorInf_Stat[n]);
return 0;
}



void S8NumericSolver::SetInitialOrder(Int_t n,Int_t a)
{
// Define criteria used to resolve ambiguity in system solution.
// Force a parameter to be larger (a=1) or smaller (a=-1) to its
// symmetric.
// n=0 -> na, nb
// n=1 -> ea1, eb1
// n=2 -> ea2, eb2
// n=3 -> ea3, eb3
// If not called, default is na > nb.
Int_t b = 1;
if(a) b=a/abs(a);
fAsym[0] = n;
fAsym[1] = b;
}

void S8NumericSolver::SetInitialOrder(std::string s1, std::string s2, std::string s3)
{
// Define criteria used to resolve ambiguity in system solution.
// Force a parameter to be larger to its
// symmetric.
// n=0 -> na, nb
// n=1 -> ea1, eb1
// n=2 -> ea2, eb2
// n=3 -> ea3, eb3
// If not called, default is na > nb.
char s = s1[1];
char t = s3[1];
char u = s1[0];
char v = s3[0];
Int_t n=-1 , a=0;
if(u=='n' && v=='n') n = 0;
if(u=='e' && v=='e')
  {
  if(s1[2]==s3[2]) n = ((Int_t) s1[2]) - 48;
  else
    {
    printf("\nIgnore SetInitialOrder(\"%s\", \"%s\", \"%s\") : Wrong parameters.\n",s1.c_str(),s2.c_str(),s3.c_str());
    return;
    }
  }
if     (s=='a' && t=='b' && s2=="<") a = -1;
else if(s=='a' && t=='b' && s2==">") a = 1;
else if(s=='b' && t=='a' && s2=="<") a = 1;
else if(s=='b' && t=='a' && s2==">") a = -1;
if(a==0 || n==-1)
  {
  printf("\nIgnore SetInitialOrder(\"%s\", \"%s\", \"%s\") : Wrong parameters.\n",s1.c_str(),s2.c_str(),s3.c_str());
  return;
  }
SetInitialOrder(n,a);
}

/********************************************************************/

Double_t S8NumericSolver::X()
{
Double_t a = V()*V()*B()*J()*U(1);
       a+= V()*V()*D()*H()*U(2);
       a+= V()*D()*J()*U(1)*U(2);
       a+= V()*V()*V()*B()*H();
     a*= kc[1][0]*kc[1][1]*kc[1][2];
       a+= kc[0][0]*kc[0][1]*kc[0][2]*(B()*H()*V()*W()*W()*-1.);
return a;
}

Double_t S8NumericSolver::Y()
{
Double_t a = U(0)*U(1)*U(2)*D()*J() ;
       a+= V()*U(1)*U(2)*(D()*I()+C()*J());
       a+= V()*U(2)*U(0)*(D()*H());
       a+= V()*U(0)*U(1)*(B()*J());
       a+= V()*V()*U(1)*(B()*I()+A()*J());
       a+= V()*V()*U(2)*(C()*H()+D()*G());
       a+= V()*V()*U(0)*(B()*H());
       a+= V()*V()*V()*(A()*H()+B()*G());
     a*= kc[1][0]*kc[1][1]*kc[1][2];
       a+= kc[0][0]*kc[0][1]*kc[0][2]*((A()*H()+B()*G())*V()*W()*W()*-1.);
       a+= -QQQ*W()*W()*D()*J();
return a;
}

Double_t S8NumericSolver::Z()
{
Double_t a = U(0)*U(1)*U(2)*(D()*I()+C()*J()) ;
       a+= V()*U(1)*U(2)*(C()*I());
       a+= V()*U(2)*U(0)*(C()*H()+D()*G());
       a+= V()*U(0)*U(1)*(B()*I()+A()*J());
       a+= V()*V()*U(1)*(A()*I());
       a+= V()*V()*U(2)*(C()*G());
       a+= V()*V()*U(0)*(A()*H()+B()*G());
       a+= V()*V()*V()*(A()*G());
     a*= kc[1][0]*kc[1][1]*kc[1][2];
       a+= kc[0][0]*kc[0][1]*kc[0][2]*(A()*G()*V()*W()*W()*-1.);
       a+= -QQQ*W()*W()*(D()*I()+C()*J());
return a;
}

Double_t S8NumericSolver::T()
{
Double_t a = U(0)*U(1)*U(2)*(C()*I()) ;
       a+= V()*U(2)*U(0)*(C()*G());
       a+= V()*U(0)*U(1)*(A()*I());
       a+= V()*V()*U(0)*(A()*G());
     a*= kc[1][0]*kc[1][1]*kc[1][2];
       a+= -QQQ*W()*W()*(C()*I());
return a;
}

Double_t S8NumericSolver::E(Int_t i,Int_t j) // i=0..1 : sample , j=0..3 : tag
{
Double_t e = -1;
if(j==0)
  {
  if(i==0) e = -V(); // first fraction
  else     e = W();
  }
else
  {
  if(i==1) e = (U(j-1)+V()*E(0,j))/W();
  else
    {
    if(j==1) e = (-K()*N()*Y()+P()*K()*X() - T()*N()*N())/((K()*N()+P()*P())*X() - P()*N()*Y() + N()*N()*Z());
    else e = (A(j-1)+B(j-1)*E(0,(j%3)+1))/(C(j-1)+D(j-1)*E(0,(j%3)+1)) ;
    }
  }
return e;
}


Double_t S8NumericSolver::ZeroFunctionFb(Double_t x)
{
fNb = x;

Double_t KK = K()*1e5;
Double_t KK2 = KK*KK;
Double_t KK3 = KK*KK2;

Double_t NN = N()*1e5;
Double_t NN2 = NN*NN;
Double_t NN3 = NN*NN2;
Double_t NN4 = NN*NN3;
Double_t NN5 = NN*NN4;

Double_t PP = P()*1e5;
Double_t PP2 = PP*PP;
Double_t PP3 = PP*PP2;

Double_t XX = X()*1e5;
Double_t XX2 = XX*XX;

Double_t YY = Y()*1e5;
Double_t YY2 = YY*YY;

Double_t ZZ = Z()*1e5;
Double_t ZZ2 = ZZ*ZZ;

Double_t TT = T()*1e5;
Double_t TT2 = TT*TT;

Double_t A1  = - KK3*NN2*XX2;
A1        += -2*KK2*NN3*XX*ZZ;
A1        += KK2*NN3*YY2   ;
A1        += -KK2*NN2*PP*XX*YY   ;
A1        += 2*KK*NN4*YY*TT   ;
A1        += -3*KK*NN3*PP*XX*TT   ;
A1        += KK*NN3*PP*YY*ZZ   ;
A1        += -KK*NN2*PP2*XX*ZZ   ;
A1        += -KK*NN4*ZZ2   ;
A1        += +NN5*TT2   ;
A1        += -NN4*PP*ZZ*TT   ;
A1        += NN3*PP2*YY*TT   ;
A1        += -NN2*PP3*XX*TT   ;

return A1;
}


Int_t S8NumericSolver::SolveSystem(Double_t* res)
{
Double_t f = 0.;
Double_t fold = 0.;
Double_t x = 0.;
Double_t xold = 0.;
Int_t n=0;

Double_t u = 1./fNpt;
for(Int_t i=0;i<=fNpt;++i)
  {
  x = i*u ;
  f = ZeroFunctionFb(x);
  if(TMath::Abs(f)<1e-6) {continue;}
  if(fold*f < 0)
    {
    if((x-xold) == 0) continue;
    Double_t aa = (f-fold)/(x-xold);
    Double_t bb = f-aa*x;
    res[n] = -bb/aa;
	// print results
    //std::cout<< "[S8Numeric] n = " << n << " res =" << res[n] << " function f = " << f << std::endl;
    n++;
    }
  fold = f;
  xold = x;
  if(n==5) break;
  }
return n;
}




