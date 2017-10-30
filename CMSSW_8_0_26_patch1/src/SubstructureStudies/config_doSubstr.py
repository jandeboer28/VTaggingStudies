####### Process initialization ##########

import FWCore.ParameterSet.Config as cms
import sys,os

process = cms.Process("Ntuple")

#Get input file name from submit script
if len(sys.argv) > 5:
	fin = sys.argv[2]
	RunEvents = int(sys.argv[3])
	JobNumber = int(sys.argv[4])
	outfolder = sys.argv[5]
	start = fin.find('MiniAODv2/')
	stop=fin.find('/MINIAODSIM')
	#fout = (os.path.split(fin)[1]).replace('.lhe','_Pythia8_'+str(JobNumber)+'_GEN.root')
	fout =  outfolder + '/' + fin[start:stop].replace("MiniAODv2","") + '_'+str(JobNumber)+'_substructure.root'
	print fout
	print '# Running input file', fin
	print '# producing outputfile', fout
	print '# Runnning job ', str(JobNumber), 'with ', str(RunEvents),' events, skipping ', str(JobNumber*RunEvents) 
	
else:
	print "No input file given"
	exit(0)
	

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_v4'

   
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(RunEvents)
)

# options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/mc/RunIISpring16MiniAODv2/WprimeToWZToWhadZhad_narrow_M-2000_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/74356808-EF40-E611-9CA9-D4AE52651C5F.root'


process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(False),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('dcap://t3se01.psi.ch:22125/'+fin)
                            )   					
							                  
####### Redo Jet clustering sequences ##########
fatjet_ptmin = 200.0

from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsPuppiSoftDrop    
from RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi import ak8PFJetsPuppiSoftDropMass
                                                                                         



#Re-running AK8 PUPPI
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.puppi.useExistingWeights = True
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')


from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets as patJetsDefault

process.load('RecoJets.JetProducers.ak8PFJetsPuppi_cfi')
process.ak8PFJetsPuppi.doAreaFastjet = True # even for standard ak8PFJets this is overwritten in RecoJets/Configuration/python/RecoPFJets_cff
from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import j2tParametersVX
process.ak8PFJetsPuppiTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
    j2tParametersVX,
    jets = cms.InputTag("ak8PFJetsPuppi")
)
process.patJetAK8PuppiCharge = cms.EDProducer("JetChargeProducer",
    src = cms.InputTag("ak8PFJetsPuppiTracksAssociatorAtVertex"),
    var = cms.string('Pt'),
    exp = cms.double(1.0)
)
    
addJetCollection(process, labelName = 'AK8Puppi',
                 jetSource = cms.InputTag('ak8PFJetsPuppi'),
                 algo= 'AK', rParam = 0.8,
				 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
				 genParticles =  cms.InputTag('prunedGenParticles'),
				 pfCandidates =  cms.InputTag('packedPFCandidates'),
                 jetCorrections = ('AK8PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 genJetCollection = cms.InputTag('slimmedGenJetsAK8')
                 )


process.load('RecoJets.JetProducers.nJettinessAdder_cfi')
process.NjettinessAK8Puppi 						= process.Njettiness.clone()
process.NjettinessAK8Puppi.src 					= cms.InputTag("ak8PFJetsPuppi")
process.NjettinessAK8Puppi.cone 				= cms.double(0.8)
process.NjettinessAK8Puppi.Njets				=cms.vuint32(1,2)
process.patJetsAK8Puppi.userData.userFloats.src = []
process.patJetsAK8Puppi.userData.userFloats.src += ['NjettinessAK8Puppi:tau1','NjettinessAK8Puppi:tau2']				 
process.selectedPatJetsAK8Puppi.cut = cms.string('pt > 200')




#Beta==0 jets --> TAGGED, used for Tau1
process.ak8PuppiJetsSoftDropBeta0 = ak8PFJetsPuppiSoftDrop.clone(jetPtMin = 150, beta = cms.double(0.0)  )
addJetCollection(process, labelName = 'AK8PuppiMMDT',
                 jetSource = cms.InputTag('ak8PuppiJetsSoftDropBeta0'),
                 algo= 'AK', rParam = 0.8,
				 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
				 genParticles =  cms.InputTag('prunedGenParticles'),
				 pfCandidates =  cms.InputTag('packedPFCandidates'),
                 jetCorrections = ('AK8PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 genJetCollection = cms.InputTag('slimmedGenJetsAK8'),
                 )
				 
				 
process.NjettinessAK8PuppiSoftdropBeta0 						= process.NjettinessAK8Puppi.clone()
process.NjettinessAK8PuppiSoftdropBeta0.src 					= cms.InputTag("ak8PuppiJetsSoftDropBeta0")
process.patJetsAK8PuppiMMDT.userData.userFloats.src = []
process.patJetsAK8PuppiMMDT.userData.userFloats.src += ['NjettinessAK8PuppiSoftdropBeta0:tau1','NjettinessAK8PuppiSoftdropBeta0:tau2']				 
process.selectedPatJetsAK8PuppiMMDT.cut = cms.string("")



#Beta==2 jets --> GROOMED, used for Tau2
process.ak8PuppiJetsSoftDropBeta2 = ak8PFJetsPuppiSoftDrop.clone(jetPtMin = 150, beta = cms.double(2.0)  )
addJetCollection(process, labelName = 'AK8PuppiSD',
                 jetSource = cms.InputTag('ak8PuppiJetsSoftDropBeta2'),
                 algo= 'AK', rParam = 0.8,
				 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
				 genParticles =  cms.InputTag('prunedGenParticles'),
				 pfCandidates =  cms.InputTag('packedPFCandidates'),
                 jetCorrections = ('AK8PFPuppi', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 genJetCollection = cms.InputTag('slimmedGenJetsAK8'),
                 )
				 
				 
process.NjettinessAK8PuppiSoftdropBeta2 						= process.NjettinessAK8Puppi.clone()
process.NjettinessAK8PuppiSoftdropBeta2.src 					= cms.InputTag("ak8PuppiJetsSoftDropBeta2")
process.patJetsAK8PuppiSD.userData.userFloats.src = []
process.patJetsAK8PuppiSD.userData.userFloats.src += ['NjettinessAK8PuppiSoftdropBeta2:tau1','NjettinessAK8PuppiSoftdropBeta2:tau2']				 
process.selectedPatJetsAK8PuppiSD.cut = cms.string("")

# Adding energy correlation functions

# For Puppi jets
from RecoJets.JetProducers.ECF_cfi import ECF 

process.ak8PFJetsPuppiECFb1 = ECF.clone(src = cms.InputTag( "ak8PFJetsPuppi" )) #Beta1
process.patJetsAK8Puppi.userData.userFloats.src += ['ak8PFJetsPuppiECFb1:ecf1','ak8PFJetsPuppiECFb1:ecf2','ak8PFJetsPuppiECFb1:ecf3']
process.ak8PFJetsPuppiECFb2 = ECF.clone( src = cms.InputTag( "ak8PFJetsPuppi" ), beta = cms.double( 2.0 ) ) #Beta2
process.patJetsAK8Puppi.userData.userFloats.src += ['ak8PFJetsPuppiECFb2:ecf1','ak8PFJetsPuppiECFb2:ecf2','ak8PFJetsPuppiECFb2:ecf3']

# For Softdrop Beta==0 jets
process.ak8PuppiJetsSoftDropBeta0ECFb1 = ECF.clone(src = cms.InputTag( "ak8PuppiJetsSoftDropBeta0" ))
process.patJetsAK8PuppiMMDT.userData.userFloats.src += ['ak8PuppiJetsSoftDropBeta0ECFb1:ecf1','ak8PuppiJetsSoftDropBeta0ECFb1:ecf2','ak8PuppiJetsSoftDropBeta0ECFb1:ecf3']
process.ak8PuppiJetsSoftDropBeta0ECFb2 = ECF.clone(src = cms.InputTag( "ak8PuppiJetsSoftDropBeta0" ), beta = cms.double( 2.0 ))
process.patJetsAK8PuppiMMDT.userData.userFloats.src += ['ak8PuppiJetsSoftDropBeta0ECFb2:ecf1','ak8PuppiJetsSoftDropBeta0ECFb2:ecf2','ak8PuppiJetsSoftDropBeta0ECFb2:ecf3']

# For Softdrop Beta==2 jets
process.ak8PuppiJetsSoftDropBeta2ECFb1 = ECF.clone(src = cms.InputTag( "ak8PuppiJetsSoftDropBeta2" ))
process.patJetsAK8PuppiSD.userData.userFloats.src += ['ak8PuppiJetsSoftDropBeta2ECFb1:ecf1','ak8PuppiJetsSoftDropBeta2ECFb1:ecf2','ak8PuppiJetsSoftDropBeta2ECFb1:ecf3']
process.ak8PuppiJetsSoftDropBeta2ECFb2 = ECF.clone(src = cms.InputTag( "ak8PuppiJetsSoftDropBeta2" ))
process.patJetsAK8PuppiSD.userData.userFloats.src += ['ak8PuppiJetsSoftDropBeta2ECFb2:ecf1','ak8PuppiJetsSoftDropBeta2ECFb2:ecf2','ak8PuppiJetsSoftDropBeta2ECFb2:ecf3']


#Assign everything to PUPPI jet
process.SoftDropBeta0ValueMap = cms.EDProducer("RecoJetToPatJetDeltaRValueMapProducer",
                                        src = cms.InputTag("ak8PFJetsPuppi"),
                                        matched = cms.InputTag("patJetsAK8PuppiMMDT"),                                         
                                        distMax = cms.double(0.8),
                                        values = cms.vstring([
                                            'userFloat("NjettinessAK8PuppiSoftdropBeta0:tau1")',
                                            'userFloat("NjettinessAK8PuppiSoftdropBeta0:tau2")',
											'userFloat("ak8PuppiJetsSoftDropBeta0ECFb1:ecf1")',
											'userFloat("ak8PuppiJetsSoftDropBeta0ECFb1:ecf2")',
											'userFloat("ak8PuppiJetsSoftDropBeta0ECFb1:ecf3")',
											'userFloat("ak8PuppiJetsSoftDropBeta0ECFb2:ecf1")',
											'userFloat("ak8PuppiJetsSoftDropBeta0ECFb2:ecf2")',
											'userFloat("ak8PuppiJetsSoftDropBeta0ECFb2:ecf3")',
                                            'pt','eta','phi','mass'
                                        ]),
                                        valueLabels = cms.vstring( [
                                            'Tau1',
                                            'Tau2',
											'ECF1b1',
											'ECF2b1',
											'ECF3b1',
											'ECF1b2',
											'ECF2b2',
											'ECF3b2',
                                            'pt','eta','phi','mass'
                                        ])
                    )
process.SoftDropBeta2ValueMap = cms.EDProducer("RecoJetToPatJetDeltaRValueMapProducer",
                                        src = cms.InputTag("ak8PFJetsPuppi"),
                                        matched = cms.InputTag("patJetsAK8PuppiSD"),                                         
                                        distMax = cms.double(0.8),
                                        values = cms.vstring([
                                            'userFloat("NjettinessAK8PuppiSoftdropBeta2:tau1")',
                                            'userFloat("NjettinessAK8PuppiSoftdropBeta2:tau2")',
											'userFloat("ak8PuppiJetsSoftDropBeta2ECFb1:ecf1")',
											'userFloat("ak8PuppiJetsSoftDropBeta2ECFb1:ecf2")',
											'userFloat("ak8PuppiJetsSoftDropBeta2ECFb1:ecf3")',
											'userFloat("ak8PuppiJetsSoftDropBeta2ECFb2:ecf1")',
											'userFloat("ak8PuppiJetsSoftDropBeta2ECFb2:ecf2")',
											'userFloat("ak8PuppiJetsSoftDropBeta2ECFb2:ecf3")',
                                            'pt','eta','phi','mass'
                                        ]),
                                        valueLabels = cms.vstring( [
                                            'Tau1',
                                            'Tau2',
											'ECF1b1',
											'ECF2b1',
											'ECF3b1',
											'ECF1b2',
											'ECF2b2',
											'ECF3b2',
                                            'pt','eta','phi','mass'
                                        ])
                    )
process.patJetsAK8Puppi.userData.userFloats.src += [cms.InputTag('SoftDropBeta0ValueMap','Tau1'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','Tau2'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','ECF1b1'),
												   cms.InputTag('SoftDropBeta0ValueMap','ECF2b1'),
												   cms.InputTag('SoftDropBeta0ValueMap','ECF3b1'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','ECF1b2'),
												   cms.InputTag('SoftDropBeta0ValueMap','ECF2b2'),
												   cms.InputTag('SoftDropBeta0ValueMap','ECF3b2'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','pt'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','eta'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','phi'),
                                                   cms.InputTag('SoftDropBeta0ValueMap','mass'),
                                                  ]
process.patJetsAK8Puppi.userData.userFloats.src += [cms.InputTag('SoftDropBeta2ValueMap','Tau1'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','Tau2'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','ECF1b1'),
												   cms.InputTag('SoftDropBeta2ValueMap','ECF2b1'),
												   cms.InputTag('SoftDropBeta2ValueMap','ECF3b1'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','ECF1b2'),
												   cms.InputTag('SoftDropBeta2ValueMap','ECF2b2'),
												   cms.InputTag('SoftDropBeta2ValueMap','ECF3b2'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','pt'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','eta'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','phi'),
                                                   cms.InputTag('SoftDropBeta2ValueMap','mass'),
                                                  ]
                                                                                



                                           
					

#Some control histograms
# process.TFileService = cms.Service("TFileService",
#     fileName = cms.string("histo.root")
# )
#
# process.plotPUPPI = cms.EDAnalyzer('CandViewHistoAnalyzer',
#     src = cms.InputTag("patJetsAK8PuppiSD"),
#     histograms = cms.VPSet(	cms.PSet(
#         						itemsToPlot = cms.untracked.int32(2),
#         						min = cms.untracked.double(0.),
#         						max = cms.untracked.double(250.),
#         						nbins = cms.untracked.int32(50),
#         						description = cms.untracked.string('jet %d PUPPI+softdrop mass [GeV]'),
#         						name = cms.untracked.string('jet_%d_mass'),
#         						plotquantity = cms.untracked.string('mass'),
# 								),
# 							)
# 						)
#
# process.plotB0 = cms.EDAnalyzer('CandViewHistoAnalyzer',
#     src = cms.InputTag("patJetsAK8PuppiMMDT"),
#     histograms = cms.VPSet(	cms.PSet(
#         						itemsToPlot = cms.untracked.int32(2),
#         						min = cms.untracked.double(0.),
#         						max = cms.untracked.double(250.),
#         						nbins = cms.untracked.int32(50),
#         						description = cms.untracked.string('jet %d PUPPI+softdrop mass [GeV]'),
#         						name = cms.untracked.string('jet_%d_mass'),
#         						plotquantity = cms.untracked.string('mass'),
# 								),
# 							)
# 						)
# process.plotB2 = cms.EDAnalyzer('CandViewHistoAnalyzer',
#     src = cms.InputTag("patJetsAK8PuppiSD"),
#     histograms = cms.VPSet(	cms.PSet(
#         						itemsToPlot = cms.untracked.int32(2),
#         						min = cms.untracked.double(0.),
#         						max = cms.untracked.double(250.),
#         						nbins = cms.untracked.int32(50),
#         						description = cms.untracked.string('jet %d PUPPI+softdrop mass [GeV]'),
#         						name = cms.untracked.string('jet_%d_mass'),
#         						plotquantity = cms.untracked.string('mass'),
# 								),
# 							)
# 						)

# process.p = cms.Path()
# process.p = process.chs*process.ak8CHSJetsSoftDrop*process.NjettinessAK8Puppi

# process.out = cms.OutputModule("PoolOutputModule",
# 		fileName=cms.untracked.string(fout),
#         outputCommands = cms.untracked.vstring('drop *',
# 					   'keep recoCandidates_*_*_*',
# 					   'keep recoCandidate_*_*_*',
# 					   'keep recoPFCandidates_*_*_*',
# 					   'keep patPackedGenParticles_*_*_*',
#    				   'keep recoGenParticles_*_*_*',
# 					   'keep recoPrunedGenParticles_*_*_*',
# 					   'keep patJets_*_*_*',
# 					   'keep recoJets_*_*_*',
# 					   'keep recoGenJets_*_*_*',
# 					   'keep *_NjettinessAK8Puppi_*_*',
# 					   'keep *_SoftDropBeta0ValueMap_*_*',
# 					   'keep *_SoftDropBeta2ValueMap_*_*',
#    					   ),
# )

# process.out = cms.OutputModule("PoolOutputModule",
# 		fileName=cms.untracked.string(fout),
#         outputCommands = cms.untracked.vstring('drop *',
# 					   'keep *_slimmedJetsAK8_*_*',
# 					   'keep *_selectedPatJetsAK8Puppi_*_*',
# 					   'keep *_selectedPatJetsAK8PuppiMMDT_*_*',
# 					   'keep *_selectedPatJetsAK8PuppiSD_*_*', 
# 					   'keep *_NjettinessAK8Puppi_*_*',
# 					   'keep *_SoftDropBeta0ValueMap_*_*',
# 					   'keep *_SoftDropBeta2ValueMap_*_*', CHECKED UNTIL HERE
# 					   'keep recoGenParticles_*_*_*',
# 					   'keep patPackedGenParticles_*_*_*',
# 					   'keep patPackedCandidates_*_*_*',
# 					   'keep recoPFJets_*_*_*', CHECKED

#
# 					  'keep *_slimmedJets_*_*',
# 					  'keep *_slimmedJetsAK8_*_*',
# 					  'keep *_slimmedJetsPuppi_*_*',
# 					  'keep *_slimmedJetsAK8PFCHSSoftDropPacked_*_*',

# 					  'keep *_slimmedGenJets_*_*',
# 					  'keep *_slimmedGenJetsAK8_*_*',

#
#
#    					   ),
# )

process.out = cms.OutputModule("PoolOutputModule",
		fileName=cms.untracked.string(fout),
        outputCommands = cms.untracked.vstring('drop *',
					   'keep patJets_slimmedJetsAK8_*_*',
					   'keep patJets_selectedPatJetsAK8Puppi_*_*',
					   # 'keep patJets_selectedPatJetsAK8PuppiMMDT_*_*',
					   # 'keep patJets_selectedPatJetsAK8PuppiSD_*_*',
					   'keep *_NjettinessAK8Puppi_*_*',
					   'keep *_SoftDropBeta0ValueMap_*_*',
					   'keep *_SoftDropBeta2ValueMap_*_*',
					   # 'keep *_ak8PFJetsPuppiValueMap_*_*',
					   # 'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_*_*',
					   # 'keep *_packedPFCandidates_*_*',
					   'keep *_puppi_*_*',   
					   'keep *_packedGenParticles_*_*',
					   'keep *_prunedGenParticles_*_*',
					   'keep *_offlineSlimmedPrimaryVertices_*_*'
					   
   					   ),
)



# process.outpath = cms.EndPath(process.out*process.plotPUPPI*process.plotB0*process.plotB2)
process.outpath = cms.EndPath(process.out)
