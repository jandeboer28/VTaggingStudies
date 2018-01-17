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
if infile.find("prime") != -1 or infile.find("Bulk") != -1:
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
jmass 		= array.array( 'f', [ 0 ] )
jsd0 		= array.array( 'f', [ 0 ] )
jsd2 		= array.array( 'f', [ 0 ] )
je 			= array.array( 'f', [ 0 ] )
jpt 		= array.array( 'f', [ 0 ] )
jeta 		= array.array( 'f', [ 0 ] )
jphi 		= array.array( 'f', [ 0 ] )
jtau21 		= array.array( 'f', [ 0 ] )
jtau21sd0 	= array.array( 'f', [ 0 ] )
jtau21sd2 	= array.array( 'f', [ 0 ] )
jN2sd0 	    = array.array( 'f', [ 0 ] )
jN2sd2 	    = array.array( 'f', [ 0 ] )
jN2b2sd0 	= array.array( 'f', [ 0 ] )
jN2b2sd2 	= array.array( 'f', [ 0 ] )
chf      	= array.array( 'f', [ 0 ] )
nhf      	= array.array( 'f', [ 0 ] )
nemf     	= array.array( 'f', [ 0 ] )
cemf     	= array.array( 'f', [ 0 ] )
chMult   	= array.array( 'f', [ 0 ] )
neMult   	= array.array( 'f', [ 0 ] )

E 			= array.array('f', n_particles*[0.])
pt 			= array.array('f', n_particles*[0.])
x 			= array.array('f', n_particles*[0.])
y 			= array.array('f', n_particles*[0.])
z 			= array.array('f', n_particles*[0.])
vtx 		= array.array('f', n_particles*[0.])
vty 		= array.array('f', n_particles*[0.])
vtz 		= array.array('f', n_particles*[0.])
nPixHits 	= array.array('f', n_particles*[0.])
fromPV 		= array.array('f', n_particles*[0.])
dxy 		= array.array('f', n_particles*[0.])
dxyErr 		= array.array('f', n_particles*[0.])
dz   		= array.array('f', n_particles*[0.])
dzErr 		= array.array('f', n_particles*[0.])
charge 		= array.array('f', n_particles*[0.])

vpdgid 		= array.array('f', [0.])
vE 			= array.array('f', [0.])
vpt 		= array.array('f', [0.])
vx 			= array.array('f', [0.])
vy 			= array.array('f', [0.])
vz 			= array.array('f', [0.])


t.Branch('vpdgid' 		, vpdgid,  		'vpdgid/f') #Real Hadronic W
t.Branch('ve' 			, vE,  			've/f') #Real Hadronic W
t.Branch('vpt' 			, vpt, 			'vpt/f')
t.Branch('vpx'			, vx,  			'vpx/f')
t.Branch('vpy'			, vy,  			'vpy/f')
t.Branch('vpz'			, vz,  			'vpz/f')
t.Branch('nconst' 		, nConst, 		'nconst/i') #Jet variables
t.Branch('jflavour' 	, jflavour, 	'jflavour/i')
t.Branch('jmass' 		, jmass, 		'jmass/f')
t.Branch('jsd0' 		, jsd0, 		'msoftdrop_beta0/f')
t.Branch('jsd2' 		, jsd2, 		'msoftdrop_beta2/f')
t.Branch('je' 			, je, 			'je/f')
t.Branch('jpt' 			, jpt, 			'jpt/f')
t.Branch('jeta' 		, jeta, 		'jeta/f')
t.Branch('jphi' 		, jphi, 		'jphi/f')
t.Branch('jchf'			, chf   ,	    'jchf/f')
t.Branch('jnhf'			, nhf   ,	    'jnhf/f')
t.Branch('jnemf'		, nemf  ,	    'jnemf/f')
t.Branch('jcemf'		, cemf  ,	    'jcemf/f')
t.Branch('jchMult'		, chMult,	    'jchMult/f')
t.Branch('jneMult'		, neMult,	    'jneMult/f')
t.Branch('jtau21' 		, jtau21, 		'jtau21/f')
t.Branch('jtau21sd0'	, jtau21sd0, 	'jtau21sd0/f')
t.Branch('jtau21sd2'	, jtau21sd2, 	'jtau21sd2/f')
t.Branch('jN2sd0'		, jN2sd0, 		'jN2sd0/f')
t.Branch('jN2sd2'		, jN2sd2, 		'jN2sd2/f')
t.Branch('jN2b2sd0'		, jN2b2sd0, 	'jN2b2sd0/f')
t.Branch('jN2b2sd2'		, jN2b2sd2, 	'jN2b2sd2/f')
t.Branch('pe' 			, E,  			'pe[%i]/f'%n_particles) #Constituent variables
t.Branch('ppt' 			, pt, 			'ppt[%i]/f'%n_particles)
t.Branch('ppx'			, x,  			'ppx[%i]/f'%n_particles)
t.Branch('ppy'			, y,  			'ppy[%i]/f'%n_particles)
t.Branch('ppz'			, z,  			'ppz[%i]/f'%n_particles)
t.Branch('ppvtx'		,   vtx 	,'ppvtx 	[%i]'%n_particles)
t.Branch('ppvty'		,   vty 	,'ppvty 	[%i]'%n_particles)
t.Branch('ppvtz'		,   vtz 	,'ppvtz 	[%i]'%n_particles)
t.Branch('ppnPixHits'   ,   nPixHits,'ppnPixHits[%i]'%n_particles)
t.Branch('ppfromPV'		,   fromPV 	,'ppfromPV 	[%i]'%n_particles)
t.Branch('ppdxy'		,   dxy 	,'ppdxy 	[%i]'%n_particles)
t.Branch('ppdxyErr'		,   dxyErr 	,'ppdxyErr 	[%i]'%n_particles)
t.Branch('ppdz'		    ,   dz   	,'ppdz   	[%i]'%n_particles)
t.Branch('ppdzErr'		,   dzErr 	,'ppdzErr 	[%i]'%n_particles)
t.Branch('ppcharge'		,   charge 	,'ppcharge 	[%i]'%n_particles)




# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

handlePruned, labelPruned  = Handle ("std::vector<reco::GenParticle>"),"prunedGenParticles"
handlePacked, labelPacked  = Handle ("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
fatjets, fatjetLabel = Handle("std::vector<pat::Jet>"), "selectedPatJetsAK8Puppi"


events = Events('dcap://t3se01.psi.ch:22125/' + infile)
# events = Events(infile)


for iev,event in enumerate(events):
	if iev % 1000 == 0:print "Event", iev
	event.getByLabel (labelPruned, handlePruned)
	pruned = handlePruned.product()
	event.getByLabel (labelPacked, handlePacked)
	packed = handlePacked.product()
	
	# Find hadronic gen Ws
	IDs = []
	genTs = []
	daughters = []
	isLep = False
	for p in pruned :
		daus = []
		if (abs(p.pdgId()) == 24 or abs(p.pdgId()) == 23)  and p.isHardProcess():
			
			for pa in pruned:
				mother = pa.mother(0)
				if mother and isAncestor(p,mother) and (abs(pa.mother(0).pdgId())==24 or abs(pa.mother(0).pdgId())==23  ):			
					if 10<abs(pa.pdgId())<19 and pa.isPromptFinalState(): 
						isLep = True
					if isLep: break
					l = ROOT.TLorentzVector() 
					l.SetPtEtaPhiM(pa.pt(),pa.eta(),pa.phi(),pa.mass())
					daus.append(l)	
			
			if isLep: continue															
			l = ROOT.TLorentzVector() 
			l.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass())
			genTs.append(l)
			IDs.append(p.pdgId())
			daughters.append(daus)
	
	if isSignal and len(genTs)==0: continue																						 
								 
	event.getByLabel(fatjetLabel, fatjets)
	nJ = 0
	for i,j in enumerate(fatjets.product()):
		nJ+= 1
		jetID = 0
		# Gen matching
		jtlv = ROOT.TLorentzVector() 
		jtlv.SetPtEtaPhiM(j.pt(),j.eta(),j.phi(),j.mass())
		isW = 0
		for gT in range(0,len(genTs)):
			dRs = 0
			for gQs in daughters[gT]:
				dR = jtlv.DeltaR(gQs)
				if dR < 0.8:
					dRs+=1
				if dRs <2:continue
				isW=1
				vpdgid[0] = IDs[gT]
				vE[0]     = genTs[gT].E()
				vpt[0]    = genTs[gT].Pt()
				vx[0]     = genTs[gT].Px()
				vy[0]     = genTs[gT].Py()
				vz[0]     = genTs[gT].Pz()
				jetID     = IDs[gT]
		
		if isSignal and not isW: continue	
		if not isSignal and isW: continue		
		
		 
		#Fill tree
		if j.userFloat("SoftDropBeta0ValueMap:Tau1") == 0.0 or j.pt()<200 or j.userFloat("SoftDropBeta0ValueMap:mass")<0.: continue
		
		jflavour 	[0] = j.pdgId()
		jtau21 		[0] = j.userFloat("NjettinessAK8Puppi:tau2")/j.userFloat("NjettinessAK8Puppi:tau1")
		jtau21sd0 	[0] = j.userFloat("SoftDropBeta0ValueMap:Tau2")/j.userFloat("SoftDropBeta0ValueMap:Tau1")
		jtau21sd2 	[0] = j.userFloat("SoftDropBeta2ValueMap:Tau2")/j.userFloat("SoftDropBeta2ValueMap:Tau1")
				
		je 				[0] = j.energy()
		jpt 			[0] = j.pt()
		jeta 			[0] = j.eta()
		jphi 			[0] = j.phi()
		jmass 			[0] = j.mass()
		jsd0            [0] = j.userFloat("SoftDropBeta0ValueMap:mass")
		jsd2            [0] = j.userFloat("SoftDropBeta2ValueMap:mass")
		
		chf    [0] = j.chargedHadronEnergyFraction()
		nhf    [0] = j.neutralHadronEnergyFraction()
		nemf   [0] = j.neutralEmEnergyFraction()
		cemf   [0] = j.chargedEmEnergyFraction()
		chMult [0] = j.chargedMultiplicity()
		neMult [0] = j.neutralMultiplicity()
		
		ecf3_norm_b1 		= j.userFloat("SoftDropBeta0ValueMap:ECF3b1")/j.userFloat("SoftDropBeta0ValueMap:ECF1b1")
		ecf2_norm_b1 		= j.userFloat("SoftDropBeta0ValueMap:ECF2b1")/j.userFloat("SoftDropBeta0ValueMap:ECF1b1")
		ecf3_norm_b2 		= j.userFloat("SoftDropBeta0ValueMap:ECF3b2")/j.userFloat("SoftDropBeta0ValueMap:ECF1b2")
		ecf2_norm_b2 		= j.userFloat("SoftDropBeta0ValueMap:ECF2b2")/j.userFloat("SoftDropBeta0ValueMap:ECF1b2")
		                    
		ecf3_norm_b1_sd     = j.userFloat("SoftDropBeta2ValueMap:ECF3b1")/j.userFloat("SoftDropBeta2ValueMap:ECF1b1")
		ecf2_norm_b1_sd     = j.userFloat("SoftDropBeta2ValueMap:ECF2b1")/j.userFloat("SoftDropBeta2ValueMap:ECF1b1")
		ecf3_norm_b2_sd     = j.userFloat("SoftDropBeta2ValueMap:ECF3b2")/j.userFloat("SoftDropBeta2ValueMap:ECF1b2")
		ecf2_norm_b2_sd     = j.userFloat("SoftDropBeta2ValueMap:ECF2b2")/j.userFloat("SoftDropBeta2ValueMap:ECF1b2")
		
		if ecf2_norm_b1==0:    jN2sd0[0] = -2
		else: 				   jN2sd0[0] = ecf3_norm_b1/(ecf2_norm_b1**2)
		if ecf2_norm_b1_sd==0: jN2sd2[0]= -2
		else: 				   jN2sd2[0] = ecf3_norm_b1_sd/(ecf2_norm_b1_sd**2)
		if ecf2_norm_b2==0:    jN2b2sd0[0]= -2
		else:                  jN2b2sd0[0] = ecf3_norm_b2/(ecf2_norm_b2**2)
		if ecf2_norm_b2_sd==0: jN2b2sd2[0]= -2
		else:                  jN2b2sd2[0] = ecf3_norm_b2_sd/(ecf2_norm_b2_sd**2)
		
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
		nConst		[0] = len(constituents)
		for i2, cand in enumerate(constituents):
			if i2>(n_particles-1): break
			
			E[i2]        = cand.energy()
			pt[i2]       = cand.pt()
			x[i2]        = cand.px()
			y[i2]        = cand.py()
			z[i2]        = cand.pz()
			vtx[i2]      = cand.vx()
			vty[i2]      = cand.vy() 
			vtz[i2]      = cand.vz() 
			nPixHits[i2] = cand.numberOfPixelHits()
			fromPV[i2]   = cand.fromPV()
			dxy[i2]      = cand.dxy()
			dxyErr[i2]   = cand.dxyError() 
			dz[i2]       = cand.dz() 
			dzErr[i2]    = cand.dzError() 
			charge[i2]   = cand.charge()

		t.Fill()
	
		

f.Write()
f.Close()   						 