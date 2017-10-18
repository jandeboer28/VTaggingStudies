# import ROOT in batch mode
import sys, os
import ROOT, array
ROOT.gROOT.SetBatch(True)

isS = sys.argv[1]

isSignal = False

if isS.find("Signal") != -1:
	f = ROOT.TFile("test_dichroic_sig.root", "recreate")
	isSignal = True
else:
	f = ROOT.TFile("test_dichroic_bkg.root", "recreate")
f.cd()
# t = ROOT.TTree("tree", "tree")

print "Is signal?  " , isSignal

class GenParticle(object):
	def __init__(self, pt, eta, phi,mass, flavour):
		self.pt = pt
		self.eta = eta
		self.phi = phi
		self.mass = mass
		self.flavour = flavour

	def get_tlv(self):
		tlv = ROOT.TLorentzVector() 
		tlv.SetPtEtaPhiM(self.pt,self.eta,self.phi,self.mass)
		return tlv
        
	def get_flavour(self):
		return self.flavour


		


def isAncestor(a,p) :
	if a == p :
		return True
	for i in xrange(0,p.numberOfMothers()) :
		if isAncestor(a,p.mother(i)) : return True
	return False
				
				
#
# jsd 		= array.array( 'f', [ 0 ] )
# je 			= array.array( 'f', [ 0 ] )
# jpt 		= array.array( 'f', [ 0 ] )
# jeta 		= array.array( 'f', [ 0 ] )
# jphi 		= array.array( 'f', [ 0 ] )
# jtau1 	= array.array( 'f', [ 0 ] )
# jtau2 	= array.array( 'f', [ 0 ] )
#
#
#
#
# t.Branch('jsd' 			, jsd, 'msoftdrop/f')
# t.Branch('je' 			, je, 'je/f')
# t.Branch('jpt' 			, jpt, 'jpt/f')
# t.Branch('jeta' 		, jeta, 'jeta/f')
# t.Branch('jphi' 		, jphi, 'jphi/f')
# t.Branch('jtau1' 		, jtau1, 'jtau1/f')
# t.Branch('jtau2' 		, jtau2, 'jtau2/f')


# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

#Handles
AK8PFJetsPuppi		, AK8PFJetsPuppiLabel	= Handle("std::vector<pat::Jet>"), "selectedPatJetsAK8Puppi"
handlePruned		, labelPruned  			= Handle ("std::vector<reco::GenParticle>"),"prunedGenParticles"
handlePacked		, labelPacked  			= Handle ("std::vector<pat::PackedGenParticle>"), "packedGenParticles"


# histsTau1 = []
# histsTau2 = []
# histsSD   = []


h_n2_sdb0	=ROOT.TH1F("h_n2_mmdt"			,""		, 100, 0., 1.) 
h_n2_sdb2	=ROOT.TH1F("h_n2_sd"			,""		, 100, 0., 1.) 
h_n2_dichr	=ROOT.TH1F("h_n2_dichr"  		,""		, 100, 0., 1.) 

h_tau21_sdb0	=ROOT.TH1F("h_tau21_mmdt"	,""		, 100, 0., 1.) 
h_tau21_sdb2	=ROOT.TH1F("h_tau21_sd"		,""		, 100, 0., 1.) 
h_tau21_dichr	=ROOT.TH1F("h_tau21_dichr"  ,""		, 100, 0., 1.) 
h_tau21_puppi	=ROOT.TH1F("h_tau21_puppi"	,""		, 100, 0., 1.) 
h_sd_sdb0		=ROOT.TH1F("h_mmdt"			,""		, 200, 0., 200) 
h_sd_sdb2		=ROOT.TH1F("h_sd"			,""		, 200, 0., 200) 
h_puppimass		=ROOT.TH1F("h_puppimass"	,""		, 200, 0., 200) 


fname = "/scratch/thaarres/SUBSTRUCTURE/Output/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp_58_substructure.root"
if isSignal: 	fname = "/scratch/thaarres/SUBSTRUCTURE/Output/WprimeToWZToWhadZhad_narrow_M-2000_13TeV-madgraph_6_substructure.root"
print "Opening file 	" ,fname
print "Saving to file 	" ,f.GetName()
events = Events(fname)
for iev,event in enumerate(events):
	if iev % 1000 == 0:print "Event", iev
	if iev > 30000: break 
	event.getByLabel (labelPruned, handlePruned)
	pruned = handlePruned.product()
	event.getByLabel (labelPacked, handlePacked)
	packed = handlePacked.product()
	
	# # Loop over gen particles, find hadronic gen Vs
	# genPs = []
	# isLep = False
	# for p in pruned :
	# 	if abs(p.pdgId()) == 23 or abs(p.pdgId()) == 24 or abs(p.pdgId()) == 25 :
	# 		for pa in packed:
	# 			mother = pa.mother(0)
	# 			if mother and isAncestor(p,mother):
	# 				if 10<abs(pa.pdgId())<19 and pa.isPromptFinalState(): isLep = True
	# 		if isLep: continue
	# 		genP = GenParticle(p.pt(),p.eta(),p.phi(),p.mass(),p.pdgId())
	# 		genPs.append(genP)
	# if isSignal and len(genPs)==0: continue
	
	
	# Loop over jets									 
	event.getByLabel(AK8PFJetsPuppiLabel, AK8PFJetsPuppi)
	for i,j in enumerate(AK8PFJetsPuppi.product()):
		if j.pt() < 200: continue
		if j.userFloat("SoftDropBeta0ValueMap:mass") < 10.0: continue
		
		# # Gen matching
# 		jtlv = ROOT.TLorentzVector()
# 		jtlv.SetPtEtaPhiM(j.pt(),j.eta(),j.phi(),j.mass())
# 		isV = 0
# 		for gV in genPs:
# 			dR = jtlv.DeltaR(gV.get_tlv())
# 			if dR < 0.8:
# 				isV = 1
# 				break
# 		if isSignal and not isV: continue
		
		puppi_tau1         = j.userFloat("NjettinessAK8Puppi:tau1")
		puppi_tau2         = j.userFloat("NjettinessAK8Puppi:tau2")
		puppi_mass		   = j.mass()
		
		beta0_tau1         = j.userFloat("SoftDropBeta0ValueMap:Tau1")
		beta0_tau2         = j.userFloat("SoftDropBeta0ValueMap:Tau2")
		beta0_sd		   = j.userFloat("SoftDropBeta0ValueMap:mass")
		
		beta2_tau1         = j.userFloat("SoftDropBeta2ValueMap:Tau1")
		beta2_tau2         = j.userFloat("SoftDropBeta2ValueMap:Tau2")
		beta2_sd		   = j.userFloat("SoftDropBeta2ValueMap:mass")
		
		puppi_ecf1		  = 10E3*j.userFloat("SoftDropBeta0ValueMap:ECF1b1")
		puppi_ecf2		  = 10E3*j.userFloat("SoftDropBeta0ValueMap:ECF2b1")
		puppi_ecf3		  = 10E3*j.userFloat("SoftDropBeta0ValueMap:ECF3b2")
		beta0_N2          = 10E3*j.userFloat("SoftDropBeta0ValueMap:ECF3b2")/(j.userFloat("SoftDropBeta0ValueMap:ECF2b1"))**2
		beta2_N2          = 10E3*j.userFloat("SoftDropBeta2ValueMap:ECF3b2")/(j.userFloat("SoftDropBeta2ValueMap:ECF2b1"))**2
		dichroic_N2		  = 10E3*(j.userFloat("SoftDropBeta2ValueMap:ECF3b2")/(j.userFloat("SoftDropBeta2ValueMap:ECF2b1"))**2)/j.userFloat("SoftDropBeta0ValueMap:ECF2b1")
		
		# print "beta0_N2     ",beta0_N2
	# 	print "beta2_N2     ",beta2_N2
	# 	print "dichroic_N2  ",dichroic_N2
	# 	print "puppi_ecf1 " ,puppi_ecf1
	# 	print "puppi_ecf2 " ,puppi_ecf2
	# 	print "puppi_ecf3 " ,puppi_ecf3
		#Softdrop jet consisting of one subjet with one constrituent, returns arbitrarly small tau1! -->Set tau21 to 1 (not two-prong like)
		
		if beta0_tau1 == 0.00000:
			print "Uh oh, tau1 is zero for MMDT tau1! tau2mmdt is = %.5f , while tau2sd is = %.5f" %(beta0_tau2,beta2_tau2)
			print "Other parameters: "
			print "mmdt Mass: " ,beta0_sd
			print "puppi Mass: " ,j.mass()
			print "sd Mass: " ,beta2_sd
			print "pt:   " ,j.userFloat("SoftDropBeta0ValueMap:pt")
			print "eta:   " ,j.userFloat("SoftDropBeta0ValueMap:eta")
			tau21_b0       = -999.
			tau21_dichroic = -999.
		else:
			tau21_b0 		= beta0_tau2/beta0_tau1
			tau21_dichroic 	= beta2_tau2/beta0_tau1
		
		tau21_b2 = beta2_tau2/beta2_tau1
		tau21_puppi = puppi_tau2/puppi_tau1
	 
		h_sd_sdb0		.Fill(beta0_sd)
		h_sd_sdb2		.Fill(beta2_sd)
		h_puppimass		.Fill(puppi_mass)
		
		if(65 <= beta0_sd <= 105.):
			h_tau21_sdb0	.Fill(tau21_b0)
			h_tau21_sdb2	.Fill(tau21_b2)
			h_tau21_dichr   .Fill(tau21_dichroic)
			h_tau21_puppi	.Fill(tau21_puppi)
			
			h_n2_sdb0	 .Fill(beta0_N2)
			h_n2_sdb0.Fill(beta2_N2)
			h_n2_dichr.Fill(dichroic_N2)

		
		
		

	# t.Fill()
	#
		

f.Write()
f.Close()   						 