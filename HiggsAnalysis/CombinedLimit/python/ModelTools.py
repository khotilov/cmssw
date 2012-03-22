import ROOT
import re, os, os.path
from sys import stderr, stdout
from math import *
ROOFIT_EXPR = "expr"
N_OBS_MAX   = 10000

class ModelBuilderBase():
    """This class defines the basic stuff for a model builder, and it's an interface on top of RooWorkspace::factory or HLF files"""
    def __init__(self,options):
        self.options = options
        self.out = stdout
        if options.bin:
            if options.out == None: options.out = re.sub(".txt$","",options.fileName)+".root"
            options.baseDir = os.path.dirname(options.fileName)
            ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
            self.out = ROOT.RooWorkspace("w","w");
            self.out._import = getattr(self.out,"import") # workaround: import is a python keyword
            self.out.dont_delete = []
            if options.verbose == 0:
                ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
            if os.environ.has_key('ROOFITSYS'):
                ROOT.gSystem.AddIncludePath(" -I%s/include " % os.environ['ROOFITSYS'])
        elif options.out != None:
            #stderr.write("Will save workspace to HLF file %s" % options.out)
            self.out = open(options.out, "w");
        if not options.bin:
            stderr.write("\nWARNING: You're not using binary mode. This is is DEPRECATED and NOT SUPPORTED anymore, and can give WRONG results.\n\n")
        if options.cexpr:
            global ROOFIT_EXPR;
            ROOFIT_EXPR = "cexpr"            
    def factory_(self,X):
        ret = self.out.factory(X);
        if ret: 
            self.out.dont_delete.append(ret)
            return ret
        else:
            print "ERROR parsing '%s'" % X
            self.out.Print("V");
            raise RuntimeError, "Error in factory statement" 
    def doComment(self,X):
        if not self.options.bin: self.out.write("// "+X+"\n");
    def doVar(self,vardef):
        if self.options.bin: self.factory_(vardef);
        else: self.out.write(vardef+";\n");
    def doSet(self,name,vars):
        if self.options.bin: self.out.defineSet(name,vars)
        else: self.out.write("%s = set(%s);\n" % (name,vars));
    def doObj(self,name,type,X):
        if self.options.bin: return self.factory_("%s::%s(%s)" % (type, name, X));
        else: self.out.write("%s = %s(%s);\n" % (name, type, X))

class ModelBuilder(ModelBuilderBase):
    """This class defines the actual methods to build a model"""
    def __init__(self,datacard,options):
        ModelBuilderBase.__init__(self,options) 
        self.DC = datacard
    def doModel(self):
        self.doObservables()
        self.doParametersOfInterest()
        self.doNuisances()
        self.doExpectedEvents()
        self.doIndividualModels()
        self.doCombination()
        if self.options.bin:
            self.doModelConfigs()
            if self.options.verbose > 1: self.out.Print("V")
            if self.options.verbose > 2: 
                print "Wrote GraphVizTree of model_s to ",self.options.out+".dot"
                self.out.pdf("model_s").graphVizTree(self.options.out+".dot", "\\n")
    def doObservables(self):
        """create pdf_bin<X> and pdf_bin<X>_bonly for each bin"""
        raise RuntimeError, "Not implemented in ModelBuilder"
    def doParametersOfInterest(self):
        self.doComment(" ----- parameters of interest -----")
        self.doComment(" --- Signal Strength --- ")
        self.doVar("r[0,20]");
        self.doComment(" --- set of all parameters of interest --- ")
        self.doSet("POI","r")
        if self.options.mass != 0:
            self.doComment(" --- We also write the higgs mass --- ")
            if self.out.var("MH"):
              self.out.var("MH").removeRange()
              self.out.var("MH").setVal(self.options.mass)
            else:
              self.doVar("MH[%g]" % self.options.mass); 
    def doNuisances(self):
        if len(self.DC.systs) == 0: return
        self.doComment(" ----- nuisances -----")
        globalobs = []
        for (n,nofloat,pdf,args,errline) in self.DC.systs: 
            if pdf == "lnN" or pdf.startswith("shape"):
                r = "-4,4" if pdf == "shape" else "-7,7"
                self.doObj("%s_Pdf" % n, "Gaussian", "%s[%s], %s_In[0,%s], 1" % (n,r,n,r));
                globalobs.append("%s_In" % n)
                if self.options.bin:
                  self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "gmM":
                val = 0;
                for c in errline.values(): #list channels
                  for v in c.values():     # list effects in each channel 
                    if v != 0:
                        if val != 0 and v != val: 
                            raise RuntimeError, "Error: line %s contains two different uncertainties %g, %g, which is not supported for gmM" % (n,v,val)
                        val = v;
                if val == 0: raise RuntimeError, "Error: line %s contains all zeroes"
                theta = val*val; kappa = 1/theta
                self.doObj("%s_Pdf" % n, "Gamma", "%s[1,%f,%f], %s_In[%g,%g,%g], %s_scaling[%g], 0" % (n, max(0.01,1-5*val), 1+5*val, n, kappa, 1, 2*kappa+4, n, theta))
                globalobs.append("%s_In" % n)
                if self.options.bin: self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "gmN":
                if False:
                    # old version, that creates a poisson with a very large range
                    self.doObj("%s_Pdf" % n, "Poisson", "%s_In[%d,0,%d], %s[0,%d], 1" % (n,args[0],2*args[0]+5,n,2*args[0]+5))
                else:
                    # new version, that creates a poisson with a narrower range (but still +/- 7 sigmas)
                    #print "Searching for bounds for",n,"poisson with obs",args[0]
                    minExp = args[0]+1 if args[0] > 0 else 0;
                    while (ROOT.TMath.Poisson(args[0], minExp) > 1e-12) and minExp > 0: 
                        #print "Poisson(%d, minExp = %f) = %g > 1e-12" % (args[0], minExp, ROOT.TMath.Poisson(args[0], minExp))
                        minExp *= 0.8;
                    maxExp = args[0]+1;
                    while (ROOT.TMath.Poisson(args[0], maxExp) > 1e-12): 
                        #print "Poisson(%d, maxExp = %f) = %g > 1e-12" % (args[0], maxExp, ROOT.TMath.Poisson(args[0], maxExp))
                        maxExp *= 1.2;
                    minObs = args[0];
                    while minObs > 0 and (ROOT.TMath.Poisson(minObs, args[0]+1) > 1e-12): 
                        #print "Poisson(minObs = %d, %f) = %g > 1e-12" % (minObs, args[0]+1, ROOT.TMath.Poisson(minObs, args[0]+1))
                        minObs -= (sqrt(args[0]) if args[0] > 10 else 1);
                    maxObs = args[0]+2;
                    while (ROOT.TMath.Poisson(maxObs, args[0]+1) > 1e-12): 
                        #print "Poisson(maxObs = %d, %f) = %g > 1e-12" % (maxObs, args[0]+1, ROOT.TMath.Poisson(maxObs, args[0]+1))
                        maxObs += (sqrt(args[0]) if args[0] > 10 else 2);
                    self.doObj("%s_Pdf" % n, "Poisson", "%s_In[%d,%f,%f], %s[%f,%f,%f], 1" % (n,args[0],minObs,maxObs,n,args[0]+1,minExp,maxExp))
                globalobs.append("%s_In" % n)
                if self.options.bin: self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "trG":
                trG_min = -7; trG_max = +7;
                for b in errline.keys():
                    for v in errline[b].values():
                        if v > 0 and 1.0 + trG_min * v < 0: trG_min = -1.0/v;
                        if v < 0 and 1.0 + trG_max * v < 0: trG_max = -1.0/v;
                r = "%f,%f" % (trG_min,trG_max);
                self.doObj("%s_Pdf" % n, "Gaussian", "%s[0,%s], %s_In[0,%s], 1" % (n,r,n,r));
                globalobs.append("%s_In" % n)
                if self.options.bin:
                  self.out.var("%s_In" % n).setConstant(True)
            elif pdf == "lnU":
                self.doObj("%s_Pdf" % n, "Uniform", "%s[-1,1]" % n);
            elif pdf == "unif":
                self.doObj("%s_Pdf" % n, "Uniform", "%s[%f,%f]" % (n,args[0],args[1]))
            elif pdf == "param":
                mean = float(args[0])
                if "/" in args[1]: 
                    sigmaL,sigmaR = args[1].split("/")
                    if sigmaL[0] != "-" or sigmaR[0] != "+": raise RuntimeError, "Asymmetric parameter uncertainties should be entered as -x/+y"
                    sigmaL = sigmaL[1:]; sigmaR = sigmaR[1:]
                    if len(args) == 3: # mean, sigma, range
                        self.doVar("%s%s" % (n,args[2]))
                    else:
                        sigma = float(args[1])
                        if self.out.var(n):
                          self.out.var(n).setConstant(False)
                          self.out.var(n).setRange(mean-4*float(sigmaL),mean+4*float(sigmaR))
                        else:
                          self.doVar("%s[%g,%g]" % (n, mean-4*float(sigmaL), mean+4*float(sigmaR)))
                    self.out.var(n).setVal(mean)
                    self.doObj("%s_Pdf" % n, "BifurGauss", "%s, %s_In[%s,%g,%g], %s, %s" % (n, n, args[0], self.out.var(n).getMin(), self.out.var(n).getMax(), sigmaL, sigmaR))
                    self.out.var("%s_In" % n).setConstant(True)
                else:
                    if len(args) == 3: # mean, sigma, range
                        self.doVar("%s%s" % (n,args[2]))
                    else:
                        sigma = float(args[1])
                        if self.out.var(n):
                          self.out.var(n).setConstant(False)
                          self.out.var(n).setRange(mean-4*sigma, mean+4*sigma)
                        else:
                          self.doVar("%s[%g,%g]" % (n, mean-4*sigma, mean+4*sigma))
                    self.out.var(n).setVal(mean)
                    self.doObj("%s_Pdf" % n, "Gaussian", "%s, %s_In[%s,%g,%g], %s" % (n, n, args[0], self.out.var(n).getMin(), self.out.var(n).getMax(), args[1]))
                    self.out.var("%s_In" % n).setConstant(True)
                globalobs.append("%s_In" % n)
            else: raise RuntimeError, "Unsupported pdf %s" % pdf
            if nofloat: 
              self.out.var("%s" % n).setAttribute("globalConstrained",True)
        if self.options.bin:
            nuisPdfs = ROOT.RooArgList()
            nuisVars = ROOT.RooArgSet()
            for (n,nf,p,a,e) in self.DC.systs:
                nuisVars.add(self.out.var(n))
                nuisPdfs.add(self.out.pdf(n+"_Pdf"))
            self.out.defineSet("nuisances", nuisVars)
            self.out.nuisPdf = ROOT.RooProdPdf("nuisancePdf", "nuisancePdf", nuisPdfs)
            self.out._import(self.out.nuisPdf)
            self.out.nuisPdfs = nuisPdfs
            gobsVars = ROOT.RooArgSet()
            for g in globalobs: gobsVars.add(self.out.var(g))
            self.out.defineSet("globalObservables", gobsVars)
        else: # doesn't work for too many nuisances :-(
            self.doSet("nuisances", ",".join(["%s"    % n for (n,nf,p,a,e) in self.DC.systs]))
            self.doObj("nuisancePdf", "PROD", ",".join(["%s_Pdf" % n for (n,nf,p,a,e) in self.DC.systs]))
            self.doSet("globalObservables", ",".join(globalobs))
    def doExpectedEvents(self):
        self.doComment(" --- Expected events in each bin, for each process ----")
        for b in self.DC.bins:
            for p in self.DC.exp[b].keys(): # so that we get only self.DC.processes contributing to this bin
                # if it's a zero background, write a zero and move on
                if self.DC.exp[b][p] == 0:
                    self.doVar("n_exp_bin%s_proc_%s[%g]" % (b, p, self.DC.exp[b][p]))
                    continue
                # collect multiplicative corrections
                strexpr = ""; strargs = ""
                gammaNorm = None; iSyst=-1
                if self.DC.isSignal[p]:
                    strexpr += " * @0";
                    strargs += ", r";
                    iSyst += 1
                for (n,nofloat,pdf,args,errline) in self.DC.systs:
                    if pdf.startswith("shape") and pdf.endswith("?"): # might be a lnN in disguise
                        if not self.isShapeSystematic(b,p,n): pdf = "lnN"
                    if pdf == "param" or pdf.startswith("shape"): continue
                    if not errline[b].has_key(p): continue
                    if errline[b][p] == 0.0: continue
                    if pdf == "lnN" and errline[b][p] == 1.0: continue
                    iSyst += 1
                    if pdf == "lnN" or pdf == "lnU":
                        if type(errline[b][p]) == list:
                            elow, ehigh = errline[b][p];
                            strexpr += " * @%d" % iSyst
                            strargs += ", AsymPow(%f,%f,%s)" % (elow, ehigh, n)
                        else:
                            strexpr += " * pow(%f,@%d)" % (errline[b][p], iSyst)
                            strargs += ", %s" % n
                    elif pdf == "gmM":
                        strexpr += " * @%d " % iSyst
                        strargs += ", %s" % n
                    elif pdf == "trG" or pdf == "unif":
                        strexpr += " * (1+%f*@%d) " % (errline[b][p], iSyst)
                        strargs += ", %s" % n
                    elif pdf == "gmN":
                        strexpr += " * @%d " % iSyst
                        strargs += ", %s" % n
                        if abs(errline[b][p] * args[0] - self.DC.exp[b][p]) > max(0.05 * max(self.DC.exp[b][p],1), errline[b][p]):
                            raise RuntimeError, "Values of N = %d, alpha = %g don't match with expected rate %g for systematics %s " % (
                                                    args[0], errline[b][p], self.DC.exp[b][p], n)
                        if gammaNorm != None:
                            raise RuntimeError, "More than one gmN uncertainty for the same bin and process (second one is %s)" % n
                        gammaNorm = "%g" % errline[b][p]
                    else: raise RuntimeError, "Unsupported pdf %s" % pdf
                # set base term (fixed or gamma)
                if gammaNorm != None:
                    strexpr = gammaNorm + strexpr
                else:
                    strexpr = str(self.DC.exp[b][p]) + strexpr
                # optimize constants
                if strargs != "":
                    self.doObj("n_exp_bin%s_proc_%s" % (b,p), ROOFIT_EXPR, "'%s'%s" % (strexpr, strargs));
                else:
                    self.doVar("n_exp_bin%s_proc_%s[%g]" % (b, p, self.DC.exp[b][p]))
    def doIndividualModels(self):
        """create pdf_bin<X> and pdf_bin<X>_bonly for each bin"""
        raise RuntimeError, "Not implemented in ModelBuilder"
    def doCombination(self):
        """create model_s and model_b pdfs"""
        raise RuntimeError, "Not implemented in ModelBuilder"
    def doModelConfigs(self):
        if not self.options.bin: raise RuntimeException
        if self.options.out == None: raise RuntimeException
        for nuis,warn in self.DC.flatParamNuisances.iteritems():
            if self.out.var(nuis): self.out.var(nuis).setAttribute("flatParam")
            elif warn: stderr.write("Missing variable %s declared as flatParam\n" % nuis)
        mc_s = ROOT.RooStats.ModelConfig("ModelConfig",       self.out)
        mc_b = ROOT.RooStats.ModelConfig("ModelConfig_bonly", self.out)
        for (l,mc) in [ ('s',mc_s), ('b',mc_b) ]:
            mc.SetPdf(self.out.pdf("model_"+l))
            mc.SetParametersOfInterest(self.out.set("POI"))
            mc.SetObservables(self.out.set("observables"))
            if len(self.DC.systs):  mc.SetNuisanceParameters(self.out.set("nuisances"))
            if self.out.set("globalObservables"): mc.SetGlobalObservables(self.out.set("globalObservables"))
            if self.options.verbose > 1: mc.Print("V")
            self.out._import(mc, mc.GetName())
        self.out.writeToFile(self.options.out)
    def isShapeSystematic(self,channel,process,syst):
        return False

class CountingModelBuilder(ModelBuilder):
    """ModelBuilder to make a counting experiment"""
    def __init__(self,datacard,options):
        ModelBuilder.__init__(self,datacard,options)
        if datacard.hasShapes: 
            raise RuntimeError, "You're using a CountingModelBuilder for a model that has shapes"
    def doObservables(self):
        if len(self.DC.obs):
            self.doComment(" ----- observables (already set to observed values) -----")
            for b in self.DC.bins: self.doVar("n_obs_bin%s[%f,0,%d]" % (b,self.DC.obs[b],N_OBS_MAX))
        else:
            self.doComment(" ----- observables -----")
            for b in self.DC.bins: self.doVar("n_obs_bin%s[0,%d]" % (b,N_OBS_MAX))
        self.doSet("observables", ",".join(["n_obs_bin%s" % b for b in self.DC.bins]))
        if len(self.DC.obs):
            if self.options.bin:
                self.out.data_obs = ROOT.RooDataSet(self.options.dataname,"observed data", self.out.set("observables"))
                self.out.data_obs.add( self.out.set("observables") )
                self.out._import(self.out.data_obs)
    def doIndividualModels(self):
        self.doComment(" --- Expected events in each bin, total (S+B and B) ----")
        for b in self.DC.bins:
            self.doObj("n_exp_bin%s_bonly" % b, "sum", ", ".join(["n_exp_bin%s_proc_%s" % (b,p) for p in self.DC.exp[b].keys() if self.DC.isSignal[p] == False]) )
            self.doObj("n_exp_bin%s"       % b, "sum", ", ".join(["n_exp_bin%s_proc_%s" % (b,p) for p in self.DC.exp[b].keys()                        ]) )
            self.doObj("pdf_bin%s"       % b, "Poisson", "n_obs_bin%s, n_exp_bin%s, 1"       % (b,b))
            self.doObj("pdf_bin%s_bonly" % b, "Poisson", "n_obs_bin%s, n_exp_bin%s_bonly, 1" % (b,b))
    def doCombination(self):
        prefix = "modelObs" if len(self.DC.systs) else "model" # if no systematics, we build directly the model
        nbins = len(self.DC.bins)
        if nbins > 50:
            from math import ceil
            nblocks = int(ceil(nbins/10.))
            for i in range(nblocks):
                self.doObj("%s_s_%d" % (prefix,i), "PROD", ",".join(["pdf_bin%s"       % self.DC.bins[j] for j in range(10*i,min(nbins,10*i+10))]))
                self.doObj("%s_b_%d" % (prefix,i), "PROD", ",".join(["pdf_bin%s_bonly" % self.DC.bins[j] for j in range(10*i,min(nbins,10*i+10))]))
            self.doObj("%s_s" % prefix, "PROD", ",".join([prefix+"_s_%d" % i for i in range(nblocks)]))
            self.doObj("%s_b" % prefix, "PROD", ",".join([prefix+"_b_%d" % i for i in range(nblocks)]))
        else: 
            self.doObj("%s_s" % prefix, "PROD", ",".join(["pdf_bin%s"       % b for b in self.DC.bins]))
            self.doObj("%s_b" % prefix, "PROD", ",".join(["pdf_bin%s_bonly" % b for b in self.DC.bins]))
        if len(self.DC.systs): # multiply by nuisances if needed
            self.doObj("model_s", "PROD", "modelObs_s, nuisancePdf")
            self.doObj("model_b", "PROD", "modelObs_b, nuisancePdf")

