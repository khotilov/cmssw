#include "Plotter.cxx"

struct Graphs {
  TGraphAsymmErrors *L1Jet_Emu;

  TGraphAsymmErrors *L1Jet10_Emu;
  TGraphAsymmErrors *L1Jet20_Emu;
  TGraphAsymmErrors *L1Jet30_Emu;
  TGraphAsymmErrors *L1Jet40_Emu;

  TGraphAsymmErrors *L1Jet10_TauVeto_Emu;
  TGraphAsymmErrors *L1Jet20_TauVeto_Emu;
  TGraphAsymmErrors *L1Jet30_TauVeto_Emu;
  TGraphAsymmErrors *L1Jet40_TauVeto_Emu;

  TGraphAsymmErrors *L1TauVeto_Emu;
  TGraphAsymmErrors *L1TauVeto_Sim;
  TGraphAsymmErrors *L1TauIsoVeto_Sim;

  void plotPF(Plotter *plotter, TString plotDir, const char *branch, const char *name, int nbins, double min, double max, TCut selection2);
  void plotL1(Plotter *plotter, TString plotDir, const char *branch, const char *name, int nbins, double min, double max);
  void combinePlots(const char *prefix, TString plotDir, const char *format, bool print);
  void fit(const char *xvar, double min, double max, TLatex& l, const char *name, TString plotdir, const char *format, bool print);
};

struct DataGraphs {
  DataGraphs(const char *filename, const char *plotdir, const char *format, bool print):
    plotter(new Plotter(filename)), plotDir(plotdir) {
    plotter->SetFormat(format);
    plotter->SetSave(print);
    this->format = format;
    this->print = print;
  }

  ~DataGraphs() {delete plotter;}
  Plotter *plotter;
  TString plotDir;
  TString format;
  bool print;

  Graphs PFTauEt;
  Graphs PFTauEta;
  Graphs PFTauEnergy;
  Graphs L1TauEt;
  Graphs L1TauEta;
  Graphs L2ClusterRMS;

  void plot();
};

void combineDataGraphs(Graphs& ztt, Graphs& qcd, TString plotDir, const char *name, const char *format, bool print);

void plotL1Efficiency_v2(bool print=false) {
  TString plotDir = "l1plots/";

  const char *format = ".C";
  //const char *format = ".png";
  //const char *format = ".gif";
  //const char *format = ".eps";
  //const char *format = ".pdf";

  DataGraphs plots_ztt("tteffAnalysis-ztautau.root", "l1plots-ztt/", format, print);
  DataGraphs plots_qcd("tteffAnalysis-qcd.root", "l1plots-qcd/", format, print);

  plots_ztt.plot();
  plots_qcd.plot();

  if(plots_ztt.PFTauEt.L1Jet_Emu && plots_qcd.PFTauEt.L1Jet_Emu) {
    combineDataGraphs(plots_ztt.PFTauEt, plots_qcd.PFTauEt, plotDir, "PFTauEt", format, print);
  }
  if(plots_ztt.PFTauEta.L1Jet_Emu && plots_qcd.PFTauEta.L1Jet_Emu) {
    combineDataGraphs(plots_ztt.PFTauEta, plots_qcd.PFTauEta, plotDir, "PFTauEta", format, print);
  }
  if(plots_ztt.PFTauEnergy.L1Jet_Emu && plots_qcd.PFTauEnergy.L1Jet_Emu) {
    combineDataGraphs(plots_ztt.PFTauEnergy, plots_qcd.PFTauEnergy, plotDir, "PFTauEnergy", format, print);
  }
  if(plots_ztt.L2ClusterRMS.L1Jet_Emu && plots_qcd.L2ClusterRMS.L1Jet_Emu) {
    combineDataGraphs(plots_ztt.L2ClusterRMS, plots_qcd.L2ClusterRMS, plotDir, "L2ClusterRMS", format, print);
  }
}

void DataGraphs::plot() {
  bool pftauet = false, pftaueta = false,  pftauenergy = false,
    l1jetet = false, l1jeteta = false, l2clusterRms = false;

  pftauet = true;
  pftaueta = true;
  //pftauenergy = true;
  l1jetet = true;
  l1jeteta = true;
  //l2clusterRms = true;

  TLatex l;
  l.SetTextSize(0.03);
  l.SetNDC(kTRUE);
  l.SetTextColor(kRed);

  gStyle->SetOptFit(1111);

  TCut DenEtaCut = "abs(PFTauEta) < 2.5";
  TCut DenEtCut = "PFTauEt>10.";

  // PFTau Et
  if(pftauet) {
    plotter->SetXTitle("PF-#tau E_{T} (GeV)");
    PFTauEt.plotPF(plotter, plotDir, "PFTauEt", "PFTauEt", 50, 0., 150., DenEtaCut);
    PFTauEt.combinePlots("PFTauEt", plotDir, format, print);
    PFTauEt.fit("E_{T}", 5., 140., l, "PFTauEt", plotDir, format, print);
  }

  // PFTau Eta
  if(pftaueta) {
    plotter->SetXTitle("PF-#tau #eta");
    PFTauEta.plotPF(plotter, plotDir, "PFTauEta", "PFTauEta", 50, -2.5, 2.5, DenEtCut);
    PFTauEta.combinePlots("PFTauEta", plotDir, format, print);
  }

  // PFTau Energy
  if(pftauenergy) {
    plotter->SetXTitle("PF-#tau energy (GeV)");
    PFTauEnergy.plotPF(plotter, plotDir, "PFTauEnergy", "PFTauEnergy", 50, 0, 200, DenEtaCut);
    PFTauEnergy.combinePlots("PFTauEnergy", plotDir, format, print);
  }

  // L2Tau Cluster RMS
  TCut l2matched = "hasMatchedL2Jet";
  if(l2clusterRms) {
    plotter->SetXTitle("L2 cluster  #DeltaR RMS");
    L2ClusterRMS.plotPF(plotter, plotDir, "L2ClusterDeltaRRMS", "L2ClusterRMS", 50, 0, 0.5, l2matched);
    L2ClusterRMS.combinePlots("L2ClusterRMS", plotDir, format, print);
  }

  // L1 Jet Et
  if(l1jetet) {
    plotter->SetXTitle("L1 jet E_{T} (GeV)");
    L1TauEt.plotL1(plotter, plotDir, "L1JetEt", "L1JetEt", 50, 0, 100);
  }

  // L1 Jet Eta
  if(l1jetet) {
    plotter->SetXTitle("L1 jet #eta");
    L1TauEta.plotL1(plotter, plotDir, "L1JetEta", "L1JetEta", 50, -2.5, 2.5);
  }
}


TLegend *combine2Plots(TGraph *plot1, TGraph *plot2,
                       const char *legend1=0, const char *legend2=0) {
  plot1->SetMinimum(0);
  plot1->SetMaximum(1.1);
  plot1->Draw("PA");
  plot1->SetMarkerColor(kBlack);
  plot2->Draw("P same");
  plot2->SetMarkerColor(kRed);

  TLegend *leg = new TLegend(0.4,0.2,0.7,0.4);
  leg->SetFillColor(kWhite);
  if(legend1 || legend2) {
    if(legend1) leg->AddEntry(plot1, legend1,"p");
    if(legend2) leg->AddEntry(plot2, legend2,"p");
    leg->Draw();
  }
  return leg;
}
TLegend *combine3Plots(TGraph *plot1, TGraph *plot2, TGraph *plot3,
                       const char *legend1=0, const char *legend2=0, const char *legend3=0) {
  TLegend *leg = combine2Plots(plot1, plot2, legend1, legend2);
  plot3->Draw("P same");
  plot3->SetMarkerColor(kBlue);

  if(legend3) {
    leg->AddEntry(plot3, legend3, "p");
    leg->Draw();
  }
  return leg;
}
TLegend *combine4Plots(TGraph *plot1, TGraph *plot2, TGraph *plot3, TGraph *plot4,
                       const char *legend1=0, const char *legend2=0, const char *legend3=0, const char *legend4=0) {
  TLegend *leg = combine3Plots(plot1, plot2, plot3, legend1, legend2, legend3);
  plot4->Draw("P same");
  plot4->SetMarkerColor(kGreen);

  if(legend4) {
    leg->AddEntry(plot4, legend4, "p");
    leg->Draw();
  }
  return leg;
}

void combineDataGraphs(Graphs& ztt, Graphs& qcd, TString plotDir, const char *name, const char *format, bool print) {
  const char *leg1 = "Z#rightarrow #tau#tau          "; // Quick&dirty: add spaces so that the text becomes smaller
  const char *leg2 = "QCD";

  combine2Plots(ztt.L1Jet_Emu, qcd.L1Jet_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet_Emu_Ztt_vs_QCD%s", name, format));

  combine2Plots(ztt.L1Jet10_Emu, qcd.L1Jet10_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet10_Emu_Ztt_vs_QCD%s", name, format));
  combine2Plots(ztt.L1Jet20_Emu, qcd.L1Jet20_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet20_Emu_Ztt_vs_QCD%s", name, format));
  combine2Plots(ztt.L1Jet30_Emu, qcd.L1Jet30_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet30_Emu_Ztt_vs_QCD%s", name, format));
  combine2Plots(ztt.L1Jet40_Emu, qcd.L1Jet40_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet40_Emu_Ztt_vs_QCD%s", name, format));

  combine2Plots(ztt.L1Jet10_TauVeto_Emu, qcd.L1Jet10_TauVeto_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1TauVeto_Jet10_Emu_Ztt_vs_QCD%s", name, format));
  combine2Plots(ztt.L1Jet20_TauVeto_Emu, qcd.L1Jet20_TauVeto_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1TauVeto_Jet20_Emu_Ztt_vs_QCD%s", name, format));
  combine2Plots(ztt.L1Jet30_TauVeto_Emu, qcd.L1Jet30_TauVeto_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1TauVeto_Jet30_Emu_Ztt_vs_QCD%s", name, format));
  combine2Plots(ztt.L1Jet40_TauVeto_Emu, qcd.L1Jet40_TauVeto_Emu, leg1, leg2);
  if(print) gPad->SaveAs(plotDir+Form("%s_L1TauVeto_Jet40_Emu_Ztt_vs_QCD%s", name, format));
}

void Graphs::combinePlots(const char *prefix, TString plotDir, const char *format, bool print) {
  // Emulator vs. Simulator
  combine2Plots(L1TauVeto_Emu, L1TauIsoVeto_Sim, "L1 Emulator", "L1 CaloSim");
  //combine2Plots(L1TauVeto_Emu, L1TauVeto_Sim, "L1 Emulator", "L1 CaloSim");
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Tau_Emu_vs_Sim%s", prefix, format));

  // Simulator veto bits
  combine3Plots(L1Jet_Emu, L1TauVeto_Sim, L1TauIsoVeto_Sim,
                "L1 Jet reconstruction",
                "L1 Shower veto",
                "L1 Shower+Isol. Veto");
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Tau_Eff%s", prefix, format));

  // Efficiency of L1 Et threshold
  combine4Plots(L1Jet10_Emu, L1Jet20_Emu, L1Jet30_Emu, L1Jet40_Emu,
                "L1 Jet E_{T} > 10",
                "L1 Jet E_{T} > 20",
                "L1 Jet E_{T} > 30",
                "L1 Jet E_{T} > 40");
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet_EtCut%s", prefix, format));

  // Efficiency of L1 tau veto
  combine4Plots(L1Jet10_TauVeto_Emu, L1Jet20_TauVeto_Emu, L1Jet30_TauVeto_Emu, L1Jet40_TauVeto_Emu, 
                "L1 Jet E_{T} > 10",
                "L1 Jet E_{T} > 20",
                "L1 Jet E_{T} > 30",
                "L1 Jet E_{T} > 40");
  if(print) gPad->SaveAs(plotDir+Form("%s_L1Tau_EtCut%s", prefix, format));
}

void Graphs::plotPF(Plotter *plotter, TString plotDir, const char *branch, const char *name, int nbins, double min, double max, TCut selection2) {
  TString draw;
  draw.Form("%s>>hnum(%d,%f,%f)", branch, nbins, min, max);

  TCut L1JetReco("hasMatchedL1Jet==1");
  TCut L1TauReco("hasMatchedL1TauJet==1");

  plotter->SetYTitle("Level-1 efficiency");
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1Jet", name));
  L1Jet_Emu = plotter->DrawHistogram(draw, L1JetReco, selection2);

  TCut L1Jet10("L1JetEt>10.");
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1Jet10", name));
  L1Jet10_Emu = plotter->DrawHistogram(draw, L1Jet10 && L1JetReco, L1JetReco && selection2);

  TCut L1Jet20("L1JetEt>20.");
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1Jet20", name));
  L1Jet20_Emu = plotter->DrawHistogram(draw, L1Jet20 && L1JetReco, L1JetReco && selection2);

  TCut L1Jet30("L1JetEt>30.");
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1Jet30", name));
  L1Jet30_Emu = plotter->DrawHistogram(draw, L1Jet30 && L1JetReco, L1JetReco && selection2);

  TCut L1Jet40("L1JetEt>40.");
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1Jet40", name));
  L1Jet40_Emu = plotter->DrawHistogram(draw, L1Jet40 && L1JetReco, L1JetReco && selection2);


  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto_Jet10", name));
  L1Jet10_TauVeto_Emu = plotter->DrawHistogram(draw, L1TauReco && L1Jet10 && L1JetReco, L1Jet10 && L1JetReco && selection2);

  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto_Jet20", name));
  L1Jet20_TauVeto_Emu = plotter->DrawHistogram(draw, L1TauReco && L1Jet20 && L1JetReco, L1Jet20 && L1JetReco && selection2);

  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto_Jet30", name));
  L1Jet30_TauVeto_Emu = plotter->DrawHistogram(draw, L1TauReco && L1Jet30 && L1JetReco, L1Jet30 && L1JetReco && selection2);

  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto_Jet40", name));
  L1Jet40_TauVeto_Emu = plotter->DrawHistogram(draw, L1TauReco && L1Jet40 && L1JetReco, L1Jet40 && L1JetReco && selection2);


  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1Tau", name));
  L1TauVeto_Emu = plotter->DrawHistogram(draw, L1TauReco, selection2);
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto", name));
  L1TauVeto_Sim = plotter->DrawHistogram(draw,"L1TauVeto==0 && hasMatchedL1Jet==1",selection2);

 
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1IsolationVeto", name));
  plotter->DrawHistogram(draw,"L1IsolationVeto==0 && hasMatchedL1Jet==1",selection2);
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto_L1IsolationVeto", name));
  L1TauIsoVeto_Sim = plotter->DrawHistogram(draw,"L1TauVeto==0&&L1IsolationVeto==0&&hasMatchedL1Jet==1",selection2);
}

void Graphs::plotL1(Plotter *plotter, TString plotDir, const char *branch, const char *name, int nbins, double min, double max) {
  TString draw;
  draw.Form("%s>>hnum(%d,%f,%f)", branch, nbins, min, max);

  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto", name));
  plotter->DrawHistogram(draw,"L1TauVeto==0");
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1IsolationVeto", name));
  plotter->DrawHistogram(draw,"L1IsolationVeto==0");
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1TauVeto_L1IsolationVeto", name));
  plotter->DrawHistogram(draw,"L1TauVeto==0&&L1IsolationVeto==0");
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1EmTauVeto", name));
  plotter->DrawHistogram(draw,"L1EmTauVeto==0");
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1HadTauVeto", name));
  plotter->DrawHistogram(draw,"L1HadTauVeto==0");
  
  plotter->SetFileName(plotDir+Form("L1Eff_%s_L1EmTauVeto_L1HadTauVeto", name));
  plotter->DrawHistogram(draw,"L1EmTauVeto==0||L1HadTauVeto==0");
}

TPaveStats *getStatsBox(TGraph *graph) {
  TList *lst = graph->GetListOfFunctions();
  if(!lst)
    return 0;

  return dynamic_cast<TPaveStats *>(lst->FindObject("stats"));
}

void fitHelperFreq(TGraph *graph, const char *xvar, double min, double max, TLatex& l) {
  graph->Draw("PA");
  graph->SetMarkerColor(kBlack);

  TF1 *myfit = new TF1("myfit","[0]*(TMath::Freq((sqrt(x)-sqrt([1]))/(2*[2])))", min, max);
  myfit->SetParameters(0,.8,1); //,.1,2,0.1);

  graph->Fit("myfit","R+");
  myfit->Draw("same");

  double xc = 0.38;
  double yc = 0.58;
  l.DrawLatex(xc, yc, Form("p_{0} #times TMath::Freq  #left(#frac{ #sqrt{%s}- #sqrt{p_{1}}}{2p_{2}}#right)", xvar)); yc -= 0.07;
  gPad->Update();
  TPaveStats *st = getStatsBox(graph);
  if(st) {
    st->SetTextSize(0.03);
    st->SetX1NDC(xc);
    st->SetY1NDC(yc-0.2);

    st->SetX2NDC(xc+0.3);
    st->SetY2NDC(yc);
  }
  else {
    l.DrawLatex(xc, yc, Form("p_{0} = %8.5f  (plateau)", myfit->GetParameter(0))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{1} = %8.5f  (turn on)", myfit->GetParameter(1))); yc -= 0.04; 
    l.DrawLatex(xc, yc, Form("p_{2} = %8.5f  (turn on)", myfit->GetParameter(2))); yc -= 0.04;
  }
}

void fitHelperFreqGaus(TGraphAsymmErrors *graph, const char *xvar, double min, double max, TLatex& l) {
  graph->Draw("PA");
  graph->SetMarkerColor(kBlack);
  TF1 *myfit = new TF1("myfit","[0]*(TMath::Freq((sqrt(x)-sqrt([1]))/(2*[2])))*TMath::Gaus(sqrt(x)-sqrt([3]),[4],[5])", min, max);
  myfit->SetParameters(0,1.,1,15.,2,1.);

  graph->Fit("myfit","R+");
  myfit->Draw("same");

  double xc = 0.3;
  double yc = 0.55;
  l.DrawLatex(xc, yc, Form("#splitline{p_{0} #times TMath::Freq  #left(#frac{ #sqrt{%s}- #sqrt{p_{1}}}{2p_{2}}#right)}{     #times TMath::Gaus  #left( #sqrt{E_{T}}- #sqrt{p_{3}},p_{4},p_{5}#right)}", xvar)); yc -= 0.07;
  gPad->Update();
  TPaveStats *st = getStatsBox(graph);
  if(st) {
    st->SetTextSize(0.03);
    st->SetX1NDC(xc);
    st->SetY1NDC(yc-0.33);

    st->SetX2NDC(xc+0.3);
    st->SetY2NDC(yc+0.02);
  }
  else {
    l.DrawLatex(xc, yc, Form("p_{0} = %8.5f  (plateau)", myfit->GetParameter(0))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{1} = %8.5f  (turn on)", myfit->GetParameter(1))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{2} = %8.5f  (turn on)", myfit->GetParameter(2))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{3} = %8.5f", myfit->GetParameter(3))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{4} = %8.5f", myfit->GetParameter(4))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{5} = %8.5f", myfit->GetParameter(5))); yc -= 0.04;
  }
}

void fitHelperLinear(TGraphAsymmErrors *graph, const char *xvar, double min, double max, TLatex& l) {
  graph->Draw("PA");
  graph->SetMarkerColor(kBlack);
  TF1 *myfit = new TF1("myfit","[0]+[1]*x", min, max);
  myfit->SetParameters(0,1.);

  graph->Fit("myfit","R+");
  myfit->Draw("same");

  double xc = 0.2;
  double yc = 0.9;
  gPad->Update();
  l.DrawLatex(xc, yc, Form("p_{0} + p_{1}%s", xvar)); yc -= 0.04;
  TPaveStats *st = getStatsBox(graph);
  if(st) {
    st->SetTextSize(0.03);
    st->SetX1NDC(xc);
    st->SetY1NDC(yc-0.2);

    st->SetX2NDC(xc+0.35);
    st->SetY2NDC(yc);
  }
  else {
    l.DrawLatex(xc, yc, Form("p_{0} = %8.5f", myfit->GetParameter(0))); yc -= 0.04;
    l.DrawLatex(xc, yc, Form("p_{1} = %8.5f", myfit->GetParameter(1))); yc -= 0.04;
  }
}

void Graphs::fit(const char *xvar, double min, double max, TLatex& l, const char *name, TString plotDir, const char *format, bool print) {
  TGraphAsymmErrors *temp = 0;

  // L1 jet efficiency
  if(L1Jet_Emu) {
    temp = dynamic_cast<TGraphAsymmErrors *>(L1Jet_Emu->Clone());
    fitHelperFreq(temp, xvar, min, max, l);
    if(print) gPad->SaveAs(plotDir+Form("%s_L1Jet_Fit%s", name, format));
  }
  /*

  // L1 tau efficiency
  if(L1TauVeto_Emu) {
    temp = dynamic_cast<TGraphAsymmErrors *>(L1TauVeto_Emu->Clone());
    fitHelperFreq(temp, xvar, min, max, l);
    if (print) gPad->SaveAs(plotDir+Form("%s_L1Tau_Fit%s",name, format));

    // Same, but with gaussian addition
    temp = dynamic_cast<TGraphAsymmErrors *>(L1TauVeto_Emu->Clone());
    fitHelperFreqGaus(temp, xvar, min, max, l);
    if (print) gPad->SaveAs(plotDir+Form("%s_L1Tau_Fit_Gaussian%s", name, format));
  }
  */
}
