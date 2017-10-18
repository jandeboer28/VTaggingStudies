# import ROOT in batch mode
import sys, os
import ROOT, array
ROOT.gROOT.SetBatch(True)

infile = sys.argv[1]
outfolder = sys.argv[2]
jobnr = sys.argv[3]

print "For job number %s :" %jobnr
print "Running over file " ,infile
print " Moving output to " , outfolder

try: os.stat(outfolder)
except: os.mkdir(outfolder)
	

isSignal = False
if infile.find("Wprime") != -1:
	isSignal  = True
	outname = "signal"
	print "Running over signal file"
else: 
	print "Running over background file"
	outname = "background"
	
f = ROOT.TFile(outfolder+"/"+outname+jobnr+".root", "recreate")
f.cd()
t = ROOT.TTree("tree", "tree")
	
print "Is signal? " , isSignal
print "Writing to outfile " , f.GetName()





def isAncestor(a,p) :
	if a == p :
		return True
	for i in xrange(0,p.numberOfMothers()) :
		if isAncestor(a,p.mother(i)) : return True
	return False
				
				
				
n_particles = 100
nConst 	= array.array( 'i', [ 0 ] )
jflavour= array.array( 'i', [ 0 ] )
jsd 		= array.array( 'f', [ 0 ] )
je 			= array.array( 'f', [ 0 ] )
jpt 		= array.array( 'f', [ 0 ] )
jeta 		= array.array( 'f', [ 0 ] )
jphi 		= array.array( 'f', [ 0 ] )
jtau1 	= array.array( 'f', [ 0 ] )
jtau2 	= array.array( 'f', [ 0 ] )
E 		= array.array('f', n_particles*[0.])
pt 		= array.array('f', n_particles*[0.])
x 		= array.array('f', n_particles*[0.])
y 		= array.array('f', n_particles*[0.])
z 		= array.array('f', n_particles*[0.])
vE 		= array.array('f', [0.])
vpt 	= array.array('f', [0.])
vx 		= array.array('f', [0.])
vy 		= array.array('f', [0.])
vz 		= array.array('f', [0.])


t.Branch('ve' 			, vE,  've/f')
t.Branch('vpt' 			, vpt, 'vpt/f')
t.Branch('vpx'			, vx,  'vpx/f')
t.Branch('vpy'			, vy,  'vpy/f')
t.Branch('vpz'			, vz,  'vpz/f')
t.Branch('nconst' 	, nConst, 'nconst/i')
t.Branch('jflavour' , jflavour, 'jflavour/i')
t.Branch('jsd' 			, jsd, 'msoftdrop/f')
t.Branch('je' 			, je, 'je/f')
t.Branch('jpt' 			, jpt, 'jpt/f')
t.Branch('jeta' 		, jeta, 'jeta/f')
t.Branch('jphi' 		, jphi, 'jphi/f')
t.Branch('jtau1' 		, jtau1, 'jtau1/f')
t.Branch('jtau2' 		, jtau2, 'jtau2/f')
t.Branch('pe' 			, E,  'pe[%i]/f'%n_particles)
t.Branch('ppt' 			, pt, 'ppt[%i]/f'%n_particles)
t.Branch('ppx'			, x,  'ppx[%i]/f'%n_particles)
t.Branch('ppy'			, y,  'ppy[%i]/f'%n_particles)
t.Branch('ppz'			, z,  'ppz[%i]/f'%n_particles)


# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

handlePruned, labelPruned  = Handle ("std::vector<reco::GenParticle>"),"prunedGenParticles"
handlePacked, labelPacked  = Handle ("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
fatjets, fatjetLabel = Handle("std::vector<pat::Jet>"), "slimmedJetsAK8"


events = Events('dcap://t3se01.psi.ch:22125/' + infile)


for iev,event in enumerate(events):
	# if iev >= 150000: break
	if iev % 1000 == 0:print "Event", iev
  
	event.getByLabel (labelPruned, handlePruned)
	pruned = handlePruned.product()
	event.getByLabel (labelPacked, handlePacked)
	packed = handlePacked.product()
	
	# Find hadronic gen Ws
	genWs = []
	isLep = False
	for p in pruned :
		if abs(p.pdgId()) == 24 :
			for pa in packed:
				mother = pa.mother(0)
				if mother and isAncestor(p,mother):
					if 10<abs(pa.pdgId())<19 and pa.isPromptFinalState(): isLep = True
			if isLep: continue															
			l = ROOT.TLorentzVector() 
			l.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass())
			genWs.append(l)
			
			
			
	if isSignal and len(genWs)==0: continue																						 
										 
	event.getByLabel(fatjetLabel, fatjets)
	nJ = 0
	for i,j in enumerate(fatjets.product()):
		nJ+= 1
		
		# Gen matching
		jtlv = ROOT.TLorentzVector() 
		jtlv.SetPtEtaPhiM(j.pt(),j.eta(),j.phi(),j.mass())
		isW = 0
		for gW in genWs:
			dR = jtlv.DeltaR(gW)
			if dR < 0.8:
				isW = 1
				vE[0] = gW.E()
				vpt[0]= gW.Pt()
				vx[0] = gW.Px()
				vy[0] = gW.Py()
				vz[0] = gW.Pz()
				break
		if isSignal and not isW: continue	
		if not isSignal and isW: continue		
		
		 
		#Fill tree
		nConst		[0] = j.numberOfDaughters()
		jflavour 	[0] = isW
		if abs(j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"))>1.: #For events with PUPPI Msd = 0.0, PUPPI TAU N = -99999
			jtau1 		[0] = -2
			jtau2 		[0] = -2
		else:
			jtau1 		[0] = j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1")
			jtau2 		[0] = j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2")
		je 				[0] = j.energy()
		jpt 			[0] = j.pt()
		jeta 			[0] = j.eta()
		jphi 			[0] = j.phi()
		
		puppi_softdrop, puppi_softdrop_subjet = ROOT.TLorentzVector() ,ROOT.TLorentzVector() 
		wSubjets = j.subjets('SoftDropPuppi')
		for iw,wsub in enumerate( wSubjets ) :
			puppi_softdrop_subjet.SetPtEtaPhiM(wsub.correctedP4(0).pt(),wsub.correctedP4(0).eta(),wsub.correctedP4(0).phi(),wsub.correctedP4(0).mass())
			puppi_softdrop+=puppi_softdrop_subjet
		jsd 			[0] = puppi_softdrop.M()

		constituents = []
		for ida in xrange( j.numberOfDaughters() ) :
			cand = j.daughter(ida)
			if cand.numberOfDaughters() == 0 :
				constituents.append( cand )
			else :
				for jda in xrange( cand.numberOfDaughters() ) :
					cand2 = cand.daughter(jda)
					constituents.append( cand2 )
		constituents.sort(key = lambda c:c.pt(), reverse=True)
		for i2, cand in enumerate(constituents):
			if i2>(n_particles-1): break
			E[i2] = cand.energy()
			pt[i2] = cand.pt()
			x[i2] = cand.px()
			y[i2] = cand.py()
			z[i2] = cand.pz()

		t.Fill()
	
		

f.Write()
f.Close()   						 