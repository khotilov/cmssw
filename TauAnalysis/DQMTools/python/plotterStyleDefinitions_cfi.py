import FWCore.ParameterSet.Config as cms
import copy

#--------------------------------------------------------------------------------
# define x-axes
#--------------------------------------------------------------------------------

xAxis = cms.PSet(
    xAxisTitle = cms.string(''),
    xAxisTitleOffset = cms.double(1.0),
    xAxisTitleSize = cms.double(0.05)
)

xAxis_pt = copy.deepcopy(xAxis)
xAxis_pt.xAxisTitle = cms.string('P_{T} / GeV')
xAxis_eta = copy.deepcopy(xAxis)
xAxis_eta.xAxisTitle = cms.string('#eta')
xAxis_phi = copy.deepcopy(xAxis)
xAxis_phi.xAxisTitle = cms.string('#phi / Rad')

xAxis_ipXY = copy.deepcopy(xAxis)
xAxis_ipXY.xAxisTitle = cms.string('IP_{xy} / mm')
xAxis_ipZ = copy.deepcopy(xAxis)
xAxis_ipZ.xAxisTitle = cms.string('IP_{z} / mm')

xAxis_dPhi = copy.deepcopy(xAxis)
xAxis_dPhi.xAxisTitle = cms.string('#Delta #phi / Rad')
xAxis_dR = copy.deepcopy(xAxis)
xAxis_dR.xAxisTitle = cms.string('#Delta R')

xAxis_prob = copy.deepcopy(xAxis)
xAxis_prob.xAxisTitle = cms.string('P')

xAxis_posX = copy.deepcopy(xAxis)
xAxis_posX.xAxisTitle = cms.string('x / cm')
xAxis_posY = copy.deepcopy(xAxis)
xAxis_posY.xAxisTitle = cms.string('y / cm')
xAxis_posR = copy.deepcopy(xAxis)
xAxis_posR.xAxisTitle = cms.string('r / cm')
xAxis_posZ = copy.deepcopy(xAxis)
xAxis_posZ.xAxisTitle = cms.string('z / cm')

xAxis_mass = copy.deepcopy(xAxis)
xAxis_mass.xAxisTitle = cms.string('M / GeV')
xAxis_transMass = copy.deepcopy(xAxis)
xAxis_transMass.xAxisTitle = cms.string('M_{T} / GeV')

xAxis_num = copy.deepcopy(xAxis)
xAxis_num.xAxisTitle = cms.string('N')

xAxis_pdgId = copy.deepcopy(xAxis)
xAxis_pdgId.xAxisTitle = cms.string('PDG Id')

xAxis_GeV = copy.deepcopy(xAxis)
xAxis_GeV.xAxisTitle = cms.string('GeV')

xAxis_unlabeled = cms.PSet(
  xAxisTitle = cms.string(''),
  xAxisTitleOffset = cms.double(1.0),
  xAxisTitleSize = cms.double(0.05)
)

#--------------------------------------------------------------------------------
# define y-axes
#--------------------------------------------------------------------------------

yAxis_numEntries_linear = cms.PSet(
    yScale = cms.string('linear'),      
    yAxisTitle = cms.string(''),        
    yAxisTitleOffset = cms.double(1.1), 
    yAxisTitleSize = cms.double(0.05)
)

yAxis_numEntries_log = cms.PSet(
    yScale = cms.string('log'), 
    minY_log = cms.double(1.e-2),
    yAxisTitle = cms.string(''), 
    yAxisTitleOffset = cms.double(1.1),
    yAxisTitleSize = cms.double(0.05)
)

yAxis_efficiency = cms.PSet(
    yScale = cms.string('linear'), 
    yAxisTitle = cms.string('#varepsilon'), 
    yAxisTitleOffset = cms.double(1.1),
    yAxisTitleSize = cms.double(0.05)
)

#--------------------------------------------------------------------------------
# define legends
#--------------------------------------------------------------------------------

legend_regular = cms.PSet(
    posX = cms.double(0.60),            
    posY = cms.double(0.64),             
    sizeX = cms.double(0.29),        
    sizeY = cms.double(0.25),            
    header = cms.string(''),          
    option = cms.string('brNDC'),       
    borderSize = cms.int32(0),          
    fillColor = cms.int32(0)             
)

legend_eff = cms.PSet(
    posX = cms.double(0.60),
    posY = cms.double(0.77),
    sizeX = cms.double(0.29),
    sizeY = cms.double(0.12),
    header = cms.string(''),
    option = cms.string('brNDC'),
    borderSize = cms.int32(0),
    fillColor = cms.int32(0)
)

#--------------------------------------------------------------------------------
# define labels
#--------------------------------------------------------------------------------

label_mcNormScale = cms.PSet(
    posX = cms.double(0.17),                
    posY = cms.double(0.77),             
    sizeX = cms.double(0.16),            
    sizeY = cms.double(0.11),          
    option = cms.string('brNDC'),   
    borderSize = cms.int32(0),         
    fillColor = cms.int32(0),           
    textColor = cms.int32(1),           
    textSize = cms.double(0.04),        
    textAlign = cms.int32(22),
    text = cms.vstring('sim. 200pb^{-1}',
                       '#sqrt{s}=10TeV')
)

#--------------------------------------------------------------------------------
# define colors
#--------------------------------------------------------------------------------

color_black = cms.int32(1)
color_red = cms.int32(628)
color_green = cms.int32(817)
color_darkBlue = cms.int32(596)
color_lightBlue = cms.int32(856)
color_yellow = cms.int32(396)
color_orange = cms.int32(797)
color_violett = cms.int32(877)
color_gray = cms.int32(15)

#--------------------------------------------------------------------------------
# define drawOptions
#--------------------------------------------------------------------------------

# for Monte Carlo expectations drawn as separate lines
# for different processes + total Standard Model expecation with error band on top
# (H1 style)
drawOption_black_separate = cms.PSet(
    lineColor = color_black,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_black,        
    fillStyle = cms.int32(0),
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_red_separate = cms.PSet(
    lineColor = color_red,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_red,        
    fillStyle = cms.int32(0),
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_green_separate = cms.PSet(
    lineColor = color_green,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_green,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_lightBlue_separate = cms.PSet(
    lineColor = color_lightBlue,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_lightBlue,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_darkBlue_separate = cms.PSet(
    lineColor = color_darkBlue,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_darkBlue,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_yellow_separate = cms.PSet(
    lineColor = color_yellow,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_yellow,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_orange_separate = cms.PSet(
    lineColor = color_orange,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_orange,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_green_separate = cms.PSet(
    lineColor = color_green,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_green,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_violett_separate = cms.PSet(
    lineColor = color_violett,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    fillColor = color_violett,
    fillStyle = cms.int32(0),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_gray_separate = cms.PSet(
    lineColor = color_gray,        
    lineStyle = cms.int32(1),         
    lineWidth = cms.int32(2),         
    fillColor = color_gray,        
    fillStyle = cms.int32(0),      
    drawOption = cms.string('hist'),  
    drawOptionLegend = cms.string('f') 
)

drawOption_sum_separate = cms.PSet(
    lineColor = color_black,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    fillColor = color_green,
    fillStyle = cms.int32(0),  
    drawOption = cms.string('eBand'),
    drawOptionLegend = cms.string('l')
)

# for Monte Carlo expectations drawn as filled areas
# stacked on top of each other
# (CDF style)      
drawOption_black_stacked = cms.PSet(
    lineColor = color_black,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_black,
    fillStyle = cms.int32(1001),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_red_stacked = cms.PSet(
    lineColor = color_red,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_red,
    fillStyle = cms.int32(3013),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_green_stacked = cms.PSet(
    lineColor = color_green,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_green,
    fillStyle = cms.int32(3004),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_darkBlue_stacked = cms.PSet(
    lineColor = color_darkBlue,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_darkBlue,
    fillStyle = cms.int32(3005),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_lightBlue_stacked = cms.PSet(
    lineColor = color_lightBlue,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_lightBlue,
    fillStyle = cms.int32(3007),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_yellow_stacked = cms.PSet(
    lineColor = color_yellow,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_yellow,
    fillStyle = cms.int32(3001),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_orange_stacked = cms.PSet(
    lineColor = color_orange,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_orange,
    fillStyle = cms.int32(3006),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_violett_stacked = cms.PSet(
    lineColor = color_violett,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_violett,
    fillStyle = cms.int32(3003),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

drawOption_gray_stacked = cms.PSet(
    lineColor = color_gray,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(1),
    fillColor = color_gray,
    fillStyle = cms.int32(1001),
    drawOption = cms.string('hist'),
    drawOptionLegend = cms.string('f')
)

# for Data (black points with error bars)
drawOption_black_points = cms.PSet(
    markerColor = color_black,
    markerSize = cms.double(1.),
    markerStyle = cms.int32(8),
    lineColor = color_black,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    drawOption = cms.string('e1p'),
    drawOptionLegend = cms.string('p')
)

# for efficiency plots
drawOption_black_eff = cms.PSet(
    markerColor = color_black,
    markerSize = cms.double(1.),
    markerStyle = cms.int32(20),
    lineColor = color_black,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    drawOption = cms.string('e1p'),
    drawOptionLegend = cms.string('p')
)

drawOption_orange_eff = cms.PSet(
    markerColor = color_orange,
    markerSize = cms.double(1.),
    markerStyle = cms.int32(20),
    lineColor = color_orange,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    drawOption = cms.string('e1p'),
    drawOptionLegend = cms.string('p')
)

drawOption_red_eff = cms.PSet(
    markerColor = color_red,
    markerSize = cms.double(1.),
    markerStyle = cms.int32(24),
    lineColor = color_red,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    drawOption = cms.string('e1p'),
    drawOptionLegend = cms.string('p')
)

drawOption_lightBlue_eff = cms.PSet(
    markerColor = color_lightBlue,
    markerSize = cms.double(1.),
    markerStyle = cms.int32(21),
    lineColor = color_lightBlue,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    drawOption = cms.string('e1p'),
    drawOptionLegend = cms.string('p')
)

drawOption_green_eff = cms.PSet(
    markerColor = color_green,
    markerSize = cms.double(1.),
    markerStyle = cms.int32(25),
    lineColor = color_green,
    lineStyle = cms.int32(1),
    lineWidth = cms.int32(2),
    drawOption = cms.string('e1p'),
    drawOptionLegend = cms.string('p')
)

#--------------------------------------------------------------------------------
# define drawOptions used per default for different Monte Carlo processes
#
# (NOTE: per default, expectations of different Monte Carlo processes
#        are stacked on top of each other)
#--------------------------------------------------------------------------------
drawOption_AH = copy.deepcopy(drawOption_yellow_stacked)
drawOption_AHbb = copy.deepcopy(drawOption_black_stacked)
drawOption_Ztautau = copy.deepcopy(drawOption_red_stacked)
drawOption_Zee = copy.deepcopy(drawOption_green_stacked)
drawOption_Zmumu = copy.deepcopy(drawOption_darkBlue_stacked)
drawOption_WplusJets = copy.deepcopy(drawOption_lightBlue_stacked)
drawOption_Wtaunu = copy.deepcopy(drawOption_red_stacked)
drawOption_Wmunu = copy.deepcopy(drawOption_darkBlue_stacked)
drawOption_Wenu = copy.deepcopy(drawOption_green_stacked)
drawOption_ZplusJets = copy.deepcopy(drawOption_lightBlue_stacked)
drawOption_ZeePlusJets = drawOption_Zee
drawOption_ZmumuPlusJets = drawOption_Zmumu
drawOption_ZtautauPlusJets = drawOption_Ztautau
drawOption_TTplusJets = copy.deepcopy(drawOption_violett_stacked)
drawOption_gammaPlusJets = copy.deepcopy(drawOption_yellow_stacked)
drawOption_QCD = copy.deepcopy(drawOption_orange_stacked)
drawOption_QCD_BCtoE_Pt20to30 = copy.deepcopy(drawOption_QCD)
drawOption_QCD_BCtoE_Pt30to80 = copy.deepcopy(drawOption_QCD)
drawOption_QCD_BCtoE_Pt80to170 = copy.deepcopy(drawOption_QCD)
drawOption_QCD_EMenriched_Pt20to30 = copy.deepcopy(drawOption_QCD)
drawOption_QCD_EMenriched_Pt30to80 = copy.deepcopy(drawOption_QCD)
drawOption_QCD_EMenriched_Pt80to170 = copy.deepcopy(drawOption_QCD)
drawOption_gammaPlusJets_Pt15to20 = copy.deepcopy(drawOption_gammaPlusJets)
drawOption_gammaPlusJets_Pt20to25 = copy.deepcopy(drawOption_gammaPlusJets)
drawOption_gammaPlusJets_Pt25to30 = copy.deepcopy(drawOption_gammaPlusJets)
drawOption_gammaPlusJets_Pt30to35 = copy.deepcopy(drawOption_gammaPlusJets)
drawOption_gammaPlusJets_PtGt35 = copy.deepcopy(drawOption_gammaPlusJets)
drawOption_gammaPlusJetsSum = copy.deepcopy(drawOption_gammaPlusJets)
drawOption_qcdSum = copy.deepcopy(drawOption_QCD)
drawOption_AH115_tautau = copy.deepcopy(drawOption_AH)
drawOption_AH160_tautau = copy.deepcopy(drawOption_AH)
drawOption_AH115bb_tautau = copy.deepcopy(drawOption_AHbb)
drawOption_AH160bb_tautau = copy.deepcopy(drawOption_AHbb)



