#
# Misc functions to manipulate Ecal records
# author: Stefano Argiro
# id: $Id: EcalCondTools.py,v 1.8 2009/07/16 08:25:49 argiro Exp $
#
#
# WARNING: we assume that the list of iovs for a given tag
#          contains one element only, the case of several elements
#          will need to be addressed

#from pluginCondDBPyInterface import *
from CondCore.Utilities import iovInspector as inspect
from ROOT import TCanvas,TH1F, TH2F, gStyle
import EcalPyUtils
import sys


def listTags(db):
    '''List all available tags for a given db '''
    tags=db.allTags()
    for tag in tags.split():
        print tag

def listIovs(db,tag):
    '''List all available iovs for a given tag'''

    try :
       iov = inspect.Iov(db,tag)
       iovlist = iov.list()
       print "Available iovs for tag: ",tag
       for p in iovlist:
           print "  Since " , p[1], " Till " , p[2]
     
    except Exception,er :
        print er 

def dumpXML(db,tag,since,filename='dump.xml'):
    '''Dump record in XML format for a given tag '''
    try :
       iov = inspect.Iov(db,tag)
       token = getToken(db,tag,since)       
       payload=inspect.PayLoad(db,token)
       out = open(filename,'w')
       print >> out, payload
      
    except Exception,er :
        print er

def plot (db, tag,since,filename='plot.root'):
    '''Invoke the plot function from the wrapper and save to the specified \
       file. The file format will reflect the extension given.'''
    
    try :
        iov = inspect.Iov(db,tag)
        iovlist = iov.list()
        token = getToken(db,tag,since)       
        payload=inspect.PayLoad(db,token)
        payload.plot(filename,"",[],[])
            
    except Exception,er :
        print er
        

def compare(tag1,db1,since1,
            tag2,db2,since2,filename='compare.root'):
  '''Produce comparison plots for two records. Save plots to file \
     according to format. tag can be an xml file'''

  coeff_1_b=[]
  coeff_2_b=[]

  coeff_1_e=[]
  coeff_2_e=[]
  
  if  tag1.find(".xml") < 0:
      try:
        exec('import '+db1.moduleName(tag1)+' as Plug')  
        what = {'how':'barrel'}
        w = inspect.setWhat(Plug.What(),what)
        ex = Plug.Extractor(w)
        p = getObject(db1,tag1,since1)
        p.extract(ex)
        coeff_1_b = [i for i in ex.values()]# first set of coefficients
        what = {'how':'endcap'}
        w = inspect.setWhat(Plug.What(),what)
        ex = Plug.Extractor(w)
        p.extract(ex)
        coeff_1_e = [i for i in ex.values()]# first set of coefficients


      except Exception,er :
          print er
  else:
      coeff_1_b,coeff_1_e = EcalPyUtils.fromXML(tag1)

  if  tag2.find(".xml")<0:
      try:  
        exec('import '+db2.moduleName(tag2)+' as Plug')
        what = {'how':'barrel'}
        w = inspect.setWhat(Plug.What(),what)
        ex = Plug.Extractor(w)
        p = getObject(db2,tag2,since2)
        p.extract(ex)
        coeff_2_b = [i for i in ex.values()]# 2nd set of coefficients
        what = {'how':'endcap'}
        w = inspect.setWhat(Plug.What(),what)
        ex = Plug.Extractor(w)
        p.extract(ex)
        coeff_2_e = [i for i in ex.values()]# first set of coefficients

        
      except Exception, er :
          print er
  else:
      coeff_2_b,coeff_2_e = EcalPyUtils.fromXML(tag2)

  gStyle.SetPalette(1)    
 

  ebhisto,ebmap, profx, profy= compareBarrel(coeff_1_b,coeff_2_b)
  eephisto,eepmap,eemhisto,eemmap=compareEndcap(coeff_1_e,coeff_2_e)

  diff_distro_can = TCanvas("difference","difference")
  diff_distro_can.Divide(2,4)

#  profx.SetMarkerType(3)
#  profy.SetMarkerType(3)
  diff_distro_can.cd(1)
  ebhisto.Draw()
  diff_distro_can.cd(2)
  ebmap.Draw("colz")
  diff_distro_can.cd(3)
  profx.Draw()
  diff_distro_can.cd(4)
  profy.Draw()
  diff_distro_can.cd(5)
  eephisto.Draw()
  diff_distro_can.cd(6)
  eepmap.Draw("colz")
  diff_distro_can.cd(7)
  eemhisto.Draw()
  diff_distro_can.cd(8)
  eemmap.Draw("colz")

  diff_distro_can.SaveAs(filename)

  eeborderphisto,eeborderpmap,eebordermhisto,eebordermmap=compareEndcapBorder(coeff_1_e,coeff_2_e)
  ebborderhisto,ebbordermap = compareBarrelBorder(coeff_1_b,coeff_2_b)

  border_diff_distro_can = TCanvas("border_difference","borders difference")
  border_diff_distro_can.Divide(2,3)

  border_diff_distro_can.cd(1)
  ebborderhisto.Draw()
  border_diff_distro_can.cd(2)
  ebbordermap.Draw("colz")
  border_diff_distro_can.cd(3)
  eeborderphisto.Draw()
  border_diff_distro_can.cd(4)
  eeborderpmap.Draw("colz")
  border_diff_distro_can.cd(5)
  eebordermhisto.Draw()
  border_diff_distro_can.cd(6)
  eebordermmap.Draw("colz")

  bordersfilename = "borders_"+filename
  prof_filename = "profiles_"+filename
  
  border_diff_distro_can.SaveAs(bordersfilename)





def histo (db, tag,since,filename='histo.root'):
    '''Make histograms and save to file. tag can be an xml file'''
    
    coeff_barl=[]
    coeff_endc=[]
    
    if  tag.find(".xml")< 0:
        try:  
          exec('import '+db.moduleName(tag)+' as Plug')

          what = {'how':'barrel'}
          w = inspect.setWhat(Plug.What(),what)
          ex = Plug.Extractor(w)
          p=getObject(db,tag,since)
          p.extract(ex)
          coeff_barl = [i for i in ex.values()]


          what = {'how':'endcap'}
          w = inspect.setWhat(Plug.What(),what)
          ex = Plug.Extractor(w)
          p.extract(ex)
          coeff_endc = [i for i in ex.values()]     

        except Exception, er :
          print er 

    else :
        coeff_barl,coeff_endc=EcalPyUtils.fromXML(tag)


    c =  TCanvas("CC distribution")
    c.Divide(2)
    eb = TH1F("EB","EB",100, -2,4)
    ee = TH1F("EE","EE",100, -2,4)

    for cb,ce in zip(coeff_barl,coeff_endc):
        eb.Fill(cb)
        ee.Fill(ce)

    c.cd(1)  
    eb.Draw()
    c.cd(2)
    ee.Draw()

    c.SaveAs(filename)


def getToken(db,tag,since):
    ''' Return payload token for a given iov, tag, db'''
    try :
       iov = inspect.Iov(db,tag)
       iovlist = iov.list()
       for p in iovlist:
           tmpsince=p[1]
           if str(tmpsince)==str(since) :
               return p[0]
       print "Could not retrieve token for tag: " , tag, " since: ", since
       sys.exit(0)
       
    except Exception, er :
       print er


def getObject(db,tag,since):
    ''' Return payload object for a given iov, tag, db'''
    found=0
    try:
       exec('import '+db.moduleName(tag)+' as Plug')  
       for elem in db.iov(tag).elements :       
           if str(elem.since())==str(since):
               found=1
               return Plug.Object(elem)
           
    except Exception, er :
        print er

    if not found :
        print "Could not retrieve payload for tag: " , tag, " since: ", since
        sys.exit(0)





def compareBarrel(coeff_barl_1,coeff_barl_2) :
  '''Return an histogram and a map of the differences '''

  diff_distro_h   = TH1F("diffh","diffh",100,-2,2)
  diff_distro_m   = TH2F("diffm","diffm",171,-85,86,360,1,361)
  diff_distro_m.SetStats(0)
  
  for i,c in enumerate(coeff_barl_1):  
      diff = c - coeff_barl_2[i]      
      ieta,iphi= EcalPyUtils.unhashEBIndex(i)
      diff_distro_h.Fill(diff) 
      diff_distro_m.Fill(ieta,iphi,diff)

  prof_x_h = diff_distro_m.ProfileX()
  prof_y_h = diff_distro_m.ProfileY()
          
  return diff_distro_h, diff_distro_m, prof_x_h, prof_y_h

def compareBarrelBorder(coeff_barl_1,coeff_barl_2) :
  '''Return an histogram and a map of the differences '''

  diff_distro_border_h   = TH1F("diffborderh","diffh",100,-2,2)
  diff_distro_border_m   = TH2F("diffborderm","diffm",171,-85,86,360,1,361)
  diff_distro_border_m.SetStats(0)
  
  for i,c in enumerate(coeff_barl_1):  
      diff = c - coeff_barl_2[i]      
      ieta,iphi= EcalPyUtils.unhashEBIndex(i)
      if (abs(ieta)==85 or abs(ieta)==65 or abs(ieta)==64 or abs(ieta)==45 or abs(ieta)==44 or abs(ieta)==25 or abs(ieta)==24 or abs(ieta)==1 or iphi%20==1 or iphi%20==0):
          diff_distro_border_h.Fill(diff) 
      if (abs(ieta)==85 or abs(ieta)==65 or abs(ieta)==64 or abs(ieta)==45 or abs(ieta)==44 or abs(ieta)==25 or abs(ieta)==24 or abs(ieta)==1 or iphi%20==0 or iphi%20==1): 
          diff_distro_border_m.Fill(ieta,iphi,diff)
          
  return diff_distro_border_h, diff_distro_border_m 




def compareEndcap(coeff_endc_1, coeff_endc_2) :
    ''' Return an histogram and a map of the differences for each endcap'''

    diff_distro_h_eep   = TH1F("diff EE+","diff EE+",100,-2,2)
    diff_distro_h_eem   = TH1F("diff EE-","diff EE-",100,-2,2)

    
    diff_distro_m_eep   = TH2F("map EE+","map EE+",100,1,101,100,1,101)
    diff_distro_m_eem   = TH2F("map EE-","map EE-",100,1,101,100,1,101)
    
    diff_distro_m_eep.SetStats(0)
    diff_distro_m_eem.SetStats(0)


    for i,c in enumerate(coeff_endc_1):  
      diff = c - coeff_endc_2[i]
      ix,iy,iz = EcalPyUtils.unhashEEIndex(i)
      
      if iz >0:
          diff_distro_h_eep.Fill(diff)
          diff_distro_m_eep.Fill(ix,iy,diff)
      else:
          diff_distro_h_eem.Fill(diff)
          diff_distro_m_eem.Fill(ix,iy,diff)

    return diff_distro_h_eep, \
           diff_distro_m_eep, \
           diff_distro_h_eem, \
           diff_distro_m_eem



def compareEndcapBorder(coeff_endc_1, coeff_endc_2) :
    ''' Return an histogram and a map of the differences for each endcap'''

    border_diff_distro_h_eep   = TH1F("borderdiff EE+","diff EE+",100,-2,2)
    border_diff_distro_h_eem   = TH1F("borderdiff EE-","diff EE-",100,-2,2)

    
    border_diff_distro_m_eep   = TH2F("bordermap EE+","map EE+",100,1,101,100,1,101)
    border_diff_distro_m_eem   = TH2F("bordermap EE-","map EE-",100,1,101,100,1,101)
    
    border_diff_distro_m_eep.SetStats(0)
    border_diff_distro_m_eem.SetStats(0)


    for i,c in enumerate(coeff_endc_1):  
      diff = c - coeff_endc_2[i]
      ix,iy,iz = EcalPyUtils.unhashEEIndex(i)
      Rsq = ((ix-50.0)**2+(iy-50.0)**2)
      
      if (iz >0 and (Rsq<144.0 or Rsq>2500.0)):
          border_diff_distro_h_eep.Fill(diff)
          border_diff_distro_m_eep.Fill(ix,iy,diff)
      elif (iz<0 and (Rsq<144.0 or Rsq>2500.0)):
          border_diff_distro_h_eem.Fill(diff)
          border_diff_distro_m_eem.Fill(ix,iy,diff)
      

    return border_diff_distro_h_eep, \
           border_diff_distro_m_eep, \
           border_diff_distro_h_eem, \
           border_diff_distro_m_eem
