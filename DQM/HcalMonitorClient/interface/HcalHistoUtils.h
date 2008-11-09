#ifndef GUARD_HCALTESTUTILS_h
#define GUARD_HCALTESTUTILS_h

#define UTILS_ETAMIN -44.5
#define UTILS_ETAMAX 44.5
#define UTILS_PHIMIN -0.5
#define UTILS_PHIMAX 73.5

// Set spacing of dashed lines in 2-D plots
#define UTILS_VERTLINESPACE 5
#define UTILS_HORIZLINESPACE 5

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

/******************************************************************************
 *
 * HcalTestUtils.h
 * v1.1
 * 4 May 2008
 * by Jeff Temple (jtemple AT fnal.gov)
 *
 *  Various template methods to grab any type of
 *  histogram stored as a MonitorElement, and 
 *  to return .gif, .html output from the histogram.
 *
 * Methods:
 *
 * myHist* getAnyHisto(myHist* hist, std::string name, std::string process, 
                         DQMStore* dbe_, bool verb, bool clone)  
 *			   
 * std::string getAnyIMG(int runNo,myHist* hist, int size, std::string htmlDir,
    		         const char* xlab, const char* ylab, bool setLogy, bool setLogx) 
 *
 * void htmlAnyHisto(int runNo, myHist *hist, 
		     const char* xlab, const char* ylab, 
		     int width, ofstream& htmlFile, 
		     std::string htmlDir, bool setLogy, bool setLogx)
 
 *
 *****************************************************************************
 */

#include "DQMServices/Core/interface/DQMStore.h"


template <class myHist>
myHist* getAnyHisto(myHist* hist,
		    std::string name, std::string process, DQMStore* dbe_,
		    bool verb, bool clone)
{
  /* 
     Method returns histogram named 'name' from DQMStore dbe_;
     'hist' pointer must be declared with 'new' (e.g., new TH2F())
     in function call so that the pointer actually points to something.
     Otherwise, the call to hist->ClassName() will fail.
  */

  using std::cout;
  using std::endl;

  if (!dbe_) return NULL;

  char title[150];  // title of histogram to grab from dbe
  char clonehisto[150];
  sprintf(title, "%sHcal/%s",process.c_str(),name.c_str());

  MonitorElement* me = dbe_->get(title); // get Monitor Element named 'title'

  if (!me) 
    {
      if (verb) cout <<"SORRY, COULD NOT FIND HISTOGRAM NAMED ["<< title<<"]"<<endl;
      return NULL; // ME not found
    } // if (!me)

  if (verb) 
    cout << "Found '" << title << "'" << endl;

  if (clone)
    sprintf(clonehisto, "ME %s",name.c_str()); // set clone histogram name


  /* As of 25 April 2008, there are 5 histogram types associated with 
     Monitor Elements (TH1F, TH2F, TH3F, TProfile, and TProfile2D).
     Provide a separate getter for each type.  Add others if necessary.
  */

  std::string histtype = hist->ClassName();

  // return TH1F from ME
  if (histtype=="TH1F")
    {
      TH1F* out;
      if (clone) out = dynamic_cast<TH1F*>(me->getTH1F()->Clone(clonehisto));
      else out = me->getTH1F();
      return dynamic_cast<myHist*>(out);
    }

  // return TH2F from ME
  else if (histtype=="TH2F")
    {
      TH2F* out;
      if (clone) out = dynamic_cast<TH2F*>(me->getTH2F()->Clone(clonehisto));
      else out = me->getTH2F();
  
      return dynamic_cast<myHist*>(out);
    }

  // return TH3F from ME
  else if (histtype=="TH3F")
    {
      TH3F* out;
      if (clone) out = dynamic_cast<TH3F*>(me->getTH3F()->Clone(clonehisto));
      else out = me->getTH3F();
      return dynamic_cast<myHist*>(out);
    }

  // return TProfile from ME
  else if (histtype=="TProfile")
    {
      TProfile* out;
      if (clone) out = dynamic_cast<TProfile*>(me->getTProfile()->Clone(clonehisto));
      else out = me->getTProfile();
      return dynamic_cast<myHist*>(out);
    }

  // return TProfile2D from ME
  else if (histtype=="TProfile2D")
    {
      TProfile2D* out;
      if (clone) out = dynamic_cast<TProfile2D*>(me->getTProfile2D()->Clone(clonehisto));
      else out = me->getTProfile2D();
      return dynamic_cast<myHist*>(out);
    }

  else
    {
      if (verb) 
	{
	  cout <<"Don't know how to access histogram '"<<title;
	  cout<<"' of type '"<<histtype<<"'"<<endl;
	}
      return NULL;
    }

  // Should never reach this point
  if (verb)
    cout <<"<HcalHistUtils::getAnyHisto>  YOU SHOULD NEVER SEE THIS MESSAGE!"<<endl;
  return NULL;

} // myHist* getAnyHisto(...)




// MAKE GIF FROM HISTOGRAM IMAGE
template <class myHist>
std::string getAnyIMG(int runNo,myHist* hist, int size, std::string htmlDir,
		      const char* xlab, const char* ylab, bool setLogy=0, bool setLogx=0 ) 
{
  /* 
     Template function draws histogram plot, and saves it as a .gif image.
     If size==1, thumbnail image is made.  Otherwise, full-size image is made.
  */

  if(hist==NULL)
    {
      return ""; // no histogram provided
    }
  
  // Grab the histogram's title, and convert it to something more palatable for use as a file name
  
  // Run cleanString algorithm  -- direct call of cleanString causes a crash 
  std::string name = (std::string)hist->GetTitle();
  //cout <<"NAME = ["<<name<<"]"<<endl;
  for ( unsigned int i = 0; i < name.size(); ++i ) {
    if ( name.substr(i, 6) == " - Run" ){
      name.replace(i, name.size()-i, "");
    }
    if ( name.substr(i, 4) == "_Run" ){
      name.replace(i, name.size()-i, "");
    }
    if ( name.substr(i, 5) == "__Run" ){
      name.replace(i, name.size()-i, "");
    }

    if (name.substr(i,1) == "(" || name.substr(i,1)==")")
      name.replace(i,1,"_");
    if (name.substr(i,1)==",")
      name.replace(i,1,"_");
    if (name.substr(i,1)=="<")
      name.replace(i,1,"_lt_");
    if (name.substr(i,1)==">")
      name.replace(i,1,"_gt_");
    if (name.substr(i,1)=="+")
      name.replace(i,1,"_plus_");
    if (name.substr(i,1)=="#")
      name.replace(i,1,"");

  } // for (unsigned int i=0; i< name.size();
  //cout <<"NEWNAME = ["<<name<<"]"<<endl;

  char dest[512]; // stores name of destination .gif file
  if(runNo>-1) sprintf(dest,"%s - Run %d",name.c_str(),runNo);
  else sprintf(dest,"%s",name.c_str());

  //hist->SetTitle(dest); // no need to change the histogram title itself, right?
  std::string title = dest;

  int xwid = 900; 
  int ywid =540;

  if(size==1) // thumbnail specified
    {
      title = title+"_tmb";
      xwid = 600; 
      ywid = 360;
    }

  // run parseString algorithm -- calling it directly causes a crash
  for ( unsigned int i = 0; i < title.size(); ++i ) {
    if ( title.substr(i, 1) == " " ){
      title.replace(i, 1, "_");
    }
    if ( title.substr(i, 1) == "#" ){
      title.replace(i, 1, "N");
    }
    if ( title.substr(i, 1) == "-" ){
      title.replace(i, 1, "_");
    }    
    if ( title.substr(i, 1) == "&" ){
      title.replace(i, 1, "_and_");
    }
    if ( title.substr(i, 1) == "(" 
	 || title.substr(i, 1) == ")" 
	 )  {
      title.replace(i, 1, "_");
    } 
    if ( title.substr(i,1) == "="){
      title.replace(i,1,"_");
    }
  } // for (unsigned int i=0; i < title.size();...)
  
  std::string outName = title+".gif";
  std::string saveName = htmlDir + outName;


  // Create canvas for histogram
  TCanvas* can = new TCanvas(dest,dest, xwid, ywid);
  hist->SetXTitle(xlab);
  hist->SetYTitle(ylab);
  std::string histtype=hist->ClassName();
 

  // Don't draw stat box for color plots
  if (((std::string)hist->GetOption())=="col" || 
      ((std::string)hist->GetOption())=="colz")
    hist->SetStats(false);

  // Draw with whatever options are set for the particular histogram

  hist->Draw(hist->GetOption());// I think that Draw should automatically use the GetOption() value, but include it here to be sure.

  // Draw Grid Lines
  int vertlinespace=UTILS_VERTLINESPACE;
  int horizlinespace=UTILS_HORIZLINESPACE;
  if (histtype=="TH2F")
    {
      // Draw vertical lines
      for (int xx=int(UTILS_ETAMIN);xx<=int(UTILS_ETAMAX);++xx)
	{
	  if (xx%vertlinespace!=0) continue;
	  TLine *vert = new TLine(xx,UTILS_PHIMIN,xx,UTILS_PHIMAX);
	  vert->SetLineStyle(3);
	  vert->Draw("same");
	}
      // Draw horizontal lines
      for (int yy=-int(UTILS_PHIMIN);yy<=int(UTILS_PHIMAX);++yy)
	{
	  if (yy%horizlinespace!=0) continue;
	  TLine *horiz = new TLine(UTILS_ETAMIN,yy,UTILS_ETAMAX,yy);
	  horiz->SetLineStyle(3);
	  horiz->Draw("same");
	}
      
    }
  
  // SetLogx, SetLogy don't seem to work.  Why not?
  if (hist->GetMaximum()>0 && hist->GetMinimum()>0)
    {
      // Don't bother with this yet until we get something useful working
      /*
      if (setLogx)
	can->SetLogx();
      if (setLogy)
	can->SetLogy();  
      */
    }	
  can->SaveAs(saveName.c_str());  
  delete can;

  return outName;
} // std::string getAnyIMG(...)




// make HTML from histogram
template <class myHist>
void htmlAnyHisto(int runNo, myHist *hist, 
		  const char* xlab, const char* ylab, 
		  int width, ofstream& htmlFile, 
		  std::string htmlDir,
		  bool setLogy=0, bool setLogx=0)
{

  /*
    Generates html output from any kind of input histogram
  */

  using std::cout;
  using std::endl;

  if(hist!=NULL)
    {    
      std::string histtype=hist->ClassName();

      // Set 2D histogram default option to "colz"
      if (histtype=="TH2F" && ((std::string)hist->GetOption())=="")
	{
	  hist->SetOption("colz");
	}

      // Form full-sized and thumbnail .gifs from histogram
      std::string imgNameTMB = "";   
      imgNameTMB = getAnyIMG(runNo,hist,1,htmlDir,xlab,ylab,setLogy, setLogx); 
      std::string imgName = "";   
      imgName = getAnyIMG(runNo,hist,2,htmlDir,xlab,ylab, setLogy, setLogx);  
      
      // Add thumbnail image to html code, linked to full-sized image
      if (imgName.size() != 0 )
	{
	htmlFile << "<td><a href=\"" <<  imgName << "\"><img src=\"" <<  imgNameTMB << "\"></a></td>" << endl;
	}
      else
	{
	  htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
	}
    } // (hist != NULL)

  else  // if no image found, make a blank table entry (maybe improve on this later? -- write an error message?)
    {
       htmlFile<<"<td align=\"center\"><br><br> Histogram does not exist in ROOT file!<br>Diagnostic flag may be off.<br>(This is normal in online running.)</td>"<<endl;
       //htmlFile << "<td><img src=\"" << " " << "\"></td>" << endl;
    }
  return;
} //void htmlAnyHisto(...)

#endif
