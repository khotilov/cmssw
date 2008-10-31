import FWCore.ParameterSet.Config as cms
from copy import deepcopy

hcalClient = cms.EDFilter("HcalMonitorClient",

                          # Variables for the Overall Client
                          runningStandalone         = cms.untracked.bool(False),
                          processName               = cms.untracked.string(''),
                          inputfile                 = cms.untracked.string(''),
                          baseHtmlDir               = cms.untracked.string('.'),
                          MonitorDaemon             = cms.untracked.bool(True),
                          diagnosticPrescaleTime    = cms.untracked.int32(-1),
                          diagnosticPrescaleEvt     = cms.untracked.int32(200),
                          diagnosticPrescaleLS      = cms.untracked.int32(-1),
                          diagnosticPrescaleUpdate  = cms.untracked.int32(-1),
                          resetFreqTime             = cms.untracked.int32(-1),
                          resetFreqEvents           = cms.untracked.int32(-1),
                          resetFreqLS               = cms.untracked.int32(-1),
                          resetFreqUpdates          = cms.untracked.int32(-1),
                          enableExit                = cms.untracked.bool(False),
                          #DoPerChanTests            = cms.untracked.bool(False), # is this used anywhere?
                          
                          # Variables from which subtasks may inherit
                          subDetsOn                 = cms.untracked.vstring('HB', 'HE', 'HF', 'HO'),
                          debug                     = cms.untracked.int32(0),
                          showTiming                = cms.untracked.bool(False),

                          # Pedestal Client,
                          PedestalClient                       = cms.untracked.bool(True),
                          PedestalClient_nominalPedMeanInADC   = cms.untracked.double(3.),
                          PedestalClient_nominalPedWidthInADC  = cms.untracked.double(1.),
                          PedestalClient_maxPedMeanDiffADC     = cms.untracked.double(1.),
                          PedestalClient_maxPedWidthDiffADC    = cms.untracked.double(1.),
                          PedestalClient_pedestalsInFC         = cms.untracked.bool(True),
                          PedestalClient_startingTimeSlice     = cms.untracked.int32(0),
                          PedestalClient_endingTimeSlice       = cms.untracked.int32(1),
                          PedestalClient_minErrorFlag          = cms.untracked.double(0.05),
                          
                          # DigiClient
                          DigiClient                = cms.untracked.bool(True),
                          digiErrorFrac             = cms.untracked.double(0.05),
                          CapIdMEAN_ErrThresh       = cms.untracked.double(1.5),
                          CapIdRMS_ErrThresh        = cms.untracked.double(0.25),

                          # Dead Cell Client
                          DeadCellClient                            = cms.untracked.bool(True),
                          DeadCellClient_test_occupancy             = cms.untracked.bool(True),
                          DeadCellClient_test_pedestal              = cms.untracked.bool(True),
                          DeadCellClient_test_energy                = cms.untracked.bool(True),
                          DeadCellClient_test_neighbor              = cms.untracked.bool(False),
                          DeadCellClient_checkNevents               = cms.untracked.int32(100),
                          DeadCellClient_checkNevents_occupancy     = cms.untracked.int32(100),
                          DeadCellClient_checkNevents_pedestal      = cms.untracked.int32(100),
                          DeadCellClient_checkNevents_energy        = cms.untracked.int32(100),
                          DeadCellClient_checkNevents_neighbor      = cms.untracked.int32(100),
                          DeadCellClient_minErrorFlag               = cms.untracked.double(0.01),
                          DeadCellClient_makeDiagnosticPlots        = cms.untracked.bool(False),

                          # Hot Cell Client
                          HotCellClient                            = cms.untracked.bool(True),
                          HotCellClient_test_persistent             = cms.untracked.bool(True),
                          HotCellClient_test_pedestal              = cms.untracked.bool(True),
                          HotCellClient_test_energy                = cms.untracked.bool(True),
                          HotCellClient_test_neighbor              = cms.untracked.bool(False),
                          HotCellClient_checkNevents               = cms.untracked.int32(100),
                          HotCellClient_checkNevents_persistent     = cms.untracked.int32(100),
                          HotCellClient_checkNevents_pedestal      = cms.untracked.int32(100),
                          HotCellClient_checkNevents_energy        = cms.untracked.int32(100),
                          HotCellClient_checkNevents_neighbor      = cms.untracked.int32(100),
                          HotCellClient_minErrorFlag               = cms.untracked.double(0.01),
                          HotCellClient_makeDiagnosticPlots        = cms.untracked.bool(False),

                          
                          
                          # DataFormatClient
                          DataFormatClient          = cms.untracked.bool(True),

                          # Summary Client
                          SummaryClient             = cms.untracked.bool(True),

                          #LED Client
                          LEDClient                 = cms.untracked.bool(True),
                          LEDRMS_ErrThresh          = cms.untracked.double(0.8),
                          LEDMEAN_ErrThresh         = cms.untracked.double(2.25),

                          # RecHit Client
                          RecHitClient              = cms.untracked.bool(True),

                          # CaloTowerClient
                          CaloTowerClient           = cms.untracked.bool(False),

                          # TrigPrimClient
                          TrigPrimClient            = cms.untracked.bool(True),

)



def setHcalClientValuesFromMonitor(client, origmonitor, debug=False):
    # need to make separate copy, or changing client after this call will also change monitor!
    monitor=deepcopy(origmonitor)
    
    #Reads variables from monitor module, and sets the client's copy of those variables to the same value.
    #This way, when you disable the DataFormat Monitor, the DataFormat client is also turned off automatically, etc.
    
    client.PedestalClient    = monitor.PedestalMonitor
    client.PedestalClient_nominalPedMeanInADC     = monitor.PedestalMonitor_nominalPedMeanInADC
    client.PedestalClient_nominalPedWidthInADC    = monitor.PedestalMonitor_nominalPedWidthInADC
    client.PedestalClient_maxPedMeanDiffADC       = monitor.PedestalMonitor_maxPedMeanDiffADC
    client.PedestalClient_maxPedWidthDiffADC      = monitor.PedestalMonitor_maxPedWidthDiffADC
    client.PedestalClient_pedestalsInFC           = monitor.PedestalMonitor_pedestalsInFC
    client.PedestalClient_startingTimeSlice       = monitor.PedestalMonitor_startingTimeSlice
    client.PedestalClient_endingTimeSlice         = monitor.PedestalMonitor_endingTimeSlice
    #client.PedestalClient_minErrorFlag            = monitor.PedestalMonitor_minErrorFlag # want to keep these separate?

    client.DeadCellClient                         = monitor.DeadCellMonitor
    client.DeadCellClient_test_occupancy          = monitor.DeadCellMonitor_test_occupancy
    client.DeadCellClient_test_pedestal           = monitor.DeadCellMonitor_test_pedestal
    client.DeadCellClient_test_energy             = monitor.DeadCellMonitor_test_energy
    client.DeadCellClient_test_neighbor           = monitor.DeadCellMonitor_test_neighbor
    client.DeadCellClient_checkNevents_occupancy  = monitor.DeadCellMonitor_checkNevents_occupancy 
    client.DeadCellClient_checkNevents_pedestal   = monitor.DeadCellMonitor_checkNevents_pedestal
    client.DeadCellClient_checkNevents_neighbor   = monitor.DeadCellMonitor_checkNevents_neighbor       
    client.DeadCellClient_checkNevents_energy     = monitor.DeadCellMonitor_checkNevents_energy        
    #client.DeadCellClient_minErrorFlag            = monitor.DeadCellMonitor_minErrorFlag # want to keep these separate?
    client.DeadCellClient_makeDiagnosticPlots     = monitor.DeadCellMonitor_makeDiagnosticPlots          


    client.HotCellClient                         = monitor.HotCellMonitor
    client.HotCellClient_test_persistent          = monitor.HotCellMonitor_test_persistent
    client.HotCellClient_test_pedestal           = monitor.HotCellMonitor_test_pedestal
    client.HotCellClient_test_energy             = monitor.HotCellMonitor_test_energy
    client.HotCellClient_test_neighbor           = monitor.HotCellMonitor_test_neighbor
    client.HotCellClient_checkNevents_persistent  = monitor.HotCellMonitor_checkNevents_persistent
    client.HotCellClient_checkNevents_pedestal   = monitor.HotCellMonitor_checkNevents_pedestal
    client.HotCellClient_checkNevents_neighbor   = monitor.HotCellMonitor_checkNevents_neighbor
    client.HotCellClient_checkNevents_energy     = monitor.HotCellMonitor_checkNevents_energy
    #client.HotCellClient_minErrorFlag            = monitor.HotCellMonitor_minErrorFlag # want to keep these separate?
    client.HotCellClient_makeDiagnosticPlots     = monitor.HotCellMonitor_makeDiagnosticPlots
                                            

    client.DigiClient        = monitor.DigiMonitor

    client.DataFormatClient  = monitor.DataFormatMonitor
    client.HotCellClient     = monitor.HotCellMonitor
    client.LEDClient         = monitor.LEDMonitor
    client.RecHitClient      = monitor.RecHitMonitor
    client.CaloTowerClient   = monitor.CaloTowerMonitor
    client.TrigPrimClient    = monitor.TrigPrimMonitor

    client.showTiming        = monitor.showTiming
    client.debug             = monitor.debug

    if (debug):
        print "HcalMonitorClient values from HcalMonitorModule: "
        print "Debug              = ", client.debug
        print "showTiming         = ", client.showTiming
        print "Pedestal Client    = ", client.PedestalClient
        print "Digi Client        = ", client.DigiClient
        print "DeadCell Client    = ", client.DeadCellClient
        print "\t\t Test DeadCell occupancy? ", client.DeadCellClient_test_occupancy
        print "\t\t Test DeadCell pedestal? ", client.DeadCellClient_test_pedestal
        print "\t\t Test DeadCell energy? ", client.DeadCellClient_test_energy
        print "\t\t Test DeadCell neighbor? ", client.DeadCellClient_test_neighbor
        print "\t\t CheckNevents DeadCell occupancy", client.DeadCellClient_checkNevents_occupancy
        print "\t\t CheckNevents DeadCell pedestal", client.DeadCellClient_checkNevents_pedestal
        print "\t\t CheckNevents DeadCell energy", client.DeadCellClient_checkNevents_energy
        print "\t\t CheckNevents DeadCell neighbor", client.DeadCellClient_checkNevents_neighbor
        print "\t\t Min Error Flag  = ",client.DeadCellClient_minErrorFlag
        print "\t\t make diagnostics? ",client.DeadCellClient_makeDiagnosticPlots

        print "HotCell Client    = ", client.HotCellClient
        print "\t\t Test HotCell persistently above threshold? ", client.HotCellClient_test_persistent
        print "\t\t Test HotCell pedestal? ", client.HotCellClient_test_pedestal
        print "\t\t Test HotCell energy? ", client.HotCellClient_test_energy
        print "\t\t Test HotCell neighbor? ", client.HotCellClient_test_neighbor
        print "\t\t CheckNevents HotCell persistent", client.HotCellClient_checkNevents_persistent
        print "\t\t CheckNevents HotCell pedestal", client.HotCellClient_checkNevents_pedestal
        print "\t\t CheckNevents HotCell energy", client.HotCellClient_checkNevents_energy
        print "\t\t CheckNevents HotCell neighbor", client.HotCellClient_checkNevents_neighbor
        print "\t\t Min Error Flag  = ",client.HotCellClient_minErrorFlag
        print "\t\t make diagnostics? ",client.HotCellClient_makeDiagnosticPlots
                                                                                        
        print "DataFormat Client  = ", client.DataFormatClient
        print "HotCell Client     = ", client.HotCellClient
        print "Summary Client     = ", client.SummaryClient
        print "LED Client         = ", client.LEDClient
        print "RecHit Client      = ", client.RecHitClient
        print "CaloTower Client   = ", client.CaloTowerClient
        print "TrigPrim Client    = ", client.TrigPrimClient

    return
