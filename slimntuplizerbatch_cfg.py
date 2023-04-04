import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')
options.register('skipEvents',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,"")
options.register('useAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,"")
options.parseArguments()

process = cms.Process('TEST')

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('input')
    )

process.maxEvents = cms.untracked.PSet(500)

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("output.root")
    )

process.ntuplizer_seq = cms.Sequence()

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
if options.useAOD : switchOnVIDElectronIdProducer(process,DataFormat.AOD)
else :           switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD)
for idmod in ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'] :
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.ntuplizer_seq = cms.Sequence()

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.load('RecoEgamma.EgammaElectronProducers.lowPtGsfElectronID_cff')
if not options.useAOD : 
    process.lowPtGsfElectronID.electrons = 'slimmedLowPtElectrons'
    process.lowPtGsfElectronID.rho = 'fixedGridRhoFastjetAll'
process.ntuplizer_seq *= process.lowPtGsfElectronID

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load('LowPtElectrons.LowPtElectrons.IDSlimNtuplizer_cfi')

process.ntuplizer_seq *= process.ntuplizer

process.ntuplizer_path = cms.Path(process.egmGsfElectronIDSequence*
                                  process.ntuplizer_seq)
process.schedule = cms.Schedule(process.ntuplizer_path)
