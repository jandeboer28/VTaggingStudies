# import ROOT in batch mode
import sys, os
import ROOT, array
ROOT.gROOT.SetBatch(True)
import numpy as np
import csv
addLola = False

n_particles_lola = 20
n_particles      = 100
n_features       = 4

if addLola:

    from modules.cola import CoLa
    from modules.lola import LoLa
    
    from keras.models import Sequential
    from keras.layers.core import Dense, Dropout, Activation, Flatten               
    from keras.models import load_model
    
    def model_lola(n_features,n_particles):

        model = Sequential() #linear stack of layers

        model.add(CoLa(input_shape = (n_features,n_particles_lola), #first layer, need to know input shape 4x20. 
                       add_total = True,
                       add_eye   = True,
                       n_out_particles = 15))

        model.add(LoLa(
            train_metric = False,
            es  = 0,
            xs  = 0,
            ys  = 0,
            zs  = 0,                 
            ms  = 1,                 
            pts = 1,                 
            n_train_es  = 1,
            n_train_ms  = 0,
            n_train_pts = 0,        
            n_train_sum_dijs   = 2,
            n_train_min_dijs   = 2))

        model.add(Flatten()) #From (4,20) to (None, 80)

        model.add(Dense(100))
        model.add(Activation('relu'))

        model.add(Dense(50))
        model.add(Activation('relu'))

        model.add(Dense(2, activation='softmax'))

        return model
    
    
    
    model = model_lola(n_features,n_particles)
    model.load_weights('/mnt/t3nfs01/data01/shome/thaarres/VTaggingStudies/CMSSW_9_4_6_patch1/src/LoLa/trainW/train_June2018/wLola_v1/wLola_v1_weights.h5')
    
    model_ptW = model_lola(n_features,n_particles)
    model.load_weights('/mnt/t3nfs01/data01/shome/thaarres/VTaggingStudies/CMSSW_9_4_6_patch1/src/LoLa/trainW/train_June2018/wLola_v1/wLola_v1_weights.h5')

def to_constit(df, n_constit, n_features):
        ret = np.expand_dims(df,axis=-1).reshape(-1, n_features, n_constit)  #expanded matrix: (1024, 4, 20)-->batch size, n_features, n_constituents
        ret = ret/500. #normalize values to be between -1 and 1 roughly
        return ret
        
image_fun = lambda x: to_constit(x,n_particles_lola, n_features)
        
infile = sys.argv[1]
outfolder = sys.argv[2]
jobnr = sys.argv[3]
ftype = sys.argv[4]
print "For job number %s :" %jobnr
print "Running over file " ,infile
print " Moving output to " , outfolder

try: os.stat(outfolder)
except: os.mkdir(outfolder)


isSignal = False
if ftype.find("S") != -1:
        isSignal  = True
        outname = "signal"
        print "Running over signal file"
else:
        print "Running over background file"
        outname = "background"

outname = infile.split("/")[8]

csvf = open(outfolder+"/"+outname+jobnr+".csv",'wb')
writer = csv.writer(csvf)

feat_list =  ["E","PX","PY","PZ"]  
brs = []
brs += ["{0}_{1}".format(feature,constit) for feature in feat_list for constit in range(n_particles)]
writer.writerow([brs])

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

type        = array.array( 'i', [ 0 ] )
i_evt       = array.array( 'i', [ 0 ] )
i_parton    = array.array( 'i', [ 0 ] )
nPV         = array.array( 'f', [ 0 ] )
nPVTrue     = array.array( 'f', [ 0 ] )
genweight   = array.array( 'f', [ 0 ] )
nConst      = array.array( 'i', [ 0 ] )
jflavour    = array.array( 'i', [ 0 ] )
jmass       = array.array( 'f', [ 0 ] )
jsd0        = array.array( 'f', [ 0 ] )
je          = array.array( 'f', [ 0 ] )
jpt         = array.array( 'f', [ 0 ] )
jeta        = array.array( 'f', [ 0 ] )
jphi        = array.array( 'f', [ 0 ] )
jtau21      = array.array( 'f', [ 0 ] )
jlola       = array.array( 'f', [ 0 ] )
jlola_ptW   = array.array( 'f', [ 0 ] )
chf         = array.array( 'f', [ 0 ] )
nhf         = array.array( 'f', [ 0 ] )
nemf        = array.array( 'f', [ 0 ] )
cemf        = array.array( 'f', [ 0 ] )
chMult      = array.array( 'f', [ 0 ] )
neMult      = array.array( 'f', [ 0 ] )

dauE        = array.array('f', 2*[0.])
daupt       = array.array('f', 2*[0.])
daueta      = array.array('f', 2*[0.])
dauphi      = array.array('f', 2*[0.])

E           = array.array('f', n_particles*[0.])
pt          = array.array('f', n_particles*[0.])
x           = array.array('f', n_particles*[0.])
y           = array.array('f', n_particles*[0.])
z           = array.array('f', n_particles*[0.])
vtx         = array.array('f', n_particles*[0.])
vty         = array.array('f', n_particles*[0.])
vtz         = array.array('f', n_particles*[0.])
nPixHits    = array.array('f', n_particles*[0.])
fromPV      = array.array('f', n_particles*[0.])
dxy         = array.array('f', n_particles*[0.])
dxyErr      = array.array('f', n_particles*[0.])
dz          = array.array('f', n_particles*[0.])
dzErr       = array.array('f', n_particles*[0.])
charge      = array.array('f', n_particles*[0.])
vpdgid      = array.array('f', [0.])
vE          = array.array('f', [0.])
vpt         = array.array('f', [0.])
veta        = array.array('f', [0.])
vphi        = array.array('f', [0.])
vx          = array.array('f', [0.])
vy          = array.array('f', [0.])
vz          = array.array('f', [0.])

t.Branch('i_parton'             , i_parton, 'i_parton/i') #number of parton in event 
t.Branch('i_evt'                , i_evt, 'i_evt/f') #Number of PV
t.Branch('type'                 , type, 'type/f') 
t.Branch('npv'                  , nPV, 'npv/f')
t.Branch('clf_lola'             , jlola, 'clf_lola/f')
t.Branch('clf_lola_ptWeighted'  , jlola, 'clf_lola_ptWeighted/f')
t.Branch('fj_corr_sdmass'       , jsd0, 'fj_corr_sdmass/f')


t.Branch('genV_pdgid'   , vpdgid,   'genV_pdgid/f') #Real Hadronic W
t.Branch('genV_e'       , vE,       'genV_e/f') #Real Hadronic W
t.Branch('genV_pt'      , vpt,      'genV_pt/f')
t.Branch('genV_eta'     , veta,     'genV_eta/f')
t.Branch('genV_phi'     , vphi,     'genV_phi/f')
t.Branch('genDaus_e'    , dauE,     'genDaus_e[2]/f') #Real Hadronic W
t.Branch('genDaus_pt'   , daupt,    'genDaus_pt[2]/f')
t.Branch('genDaus_eta'  , daueta,   'genDaus_eta[2]/f')
t.Branch('genDaus_phi'  , dauphi,   'genDaus_phi[2]/f')

t.Branch('genWeight'    , genweight,'genWeight/f') #Real Hadronic W
t.Branch('nconst'       , nConst,   'nconst/i') #Jet variables
t.Branch('jmass'        , jmass,    'jmass/f')

t.Branch('je'           , je,       'je/f')
t.Branch('jpt'          , jpt,      'jpt/f')
t.Branch('jeta'         , jeta,     'jeta/f')
t.Branch('jphi'         , jphi,     'jphi/f')
t.Branch('jchf'         , chf,      'jchf/f')
t.Branch('jnhf'         , nhf,      'jnhf/f')
t.Branch('jnemf'        , nemf,     'jnemf/f')
t.Branch('jcemf'        , cemf,     'jcemf/f')
t.Branch('jchMult'      , chMult,   'jchMult/f')
t.Branch('jneMult'      , neMult,   'jneMult/f')
t.Branch('jtau21'       , jtau21,   'jtau21/f')

t.Branch('pe'           , E,        'pe[%i]/f'%n_particles) #Constituent variables
t.Branch('ppt'          , pt,       'ppt[%i]/f'%n_particles)
t.Branch('ppx'          , x,        'ppx[%i]/f'%n_particles)
t.Branch('ppy'          , y,        'ppy[%i]/f'%n_particles)
t.Branch('ppz'          , z,        'ppz[%i]/f'%n_particles)
t.Branch('ppvtx'        , vtx,      'ppvtx[%i]'%n_particles)
t.Branch('ppvty'        , vty,      'ppvty[%i]'%n_particles)
t.Branch('ppvtz'        , vtz,      'ppvtz[%i]'%n_particles)
t.Branch('ppnPixHits'   , nPixHits, 'ppnPixHits[%i]'%n_particles)
t.Branch('ppfromPV'     , fromPV,   'ppfromPV[%i]'%n_particles)
t.Branch('ppdxy'        , dxy,      'ppdxy[%i]'%n_particles)
t.Branch('ppdxyErr'     , dxyErr,   'ppdxyErr[%i]'%n_particles)
t.Branch('ppdz'         , dz,       'ppdz[%i]'%n_particles)
t.Branch('ppdzErr'      , dzErr,    'ppdzErr[%i]'%n_particles)
t.Branch('ppcharge'     , charge,   'ppcharge[%i]'%n_particles)




# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
pvInfo, pvLabel      = Handle(" std::vector<PileupSummaryInfo>"), "slimmedAddPileupInfo"
genInfo, genInfoLabel= Handle("GenEventInfoProduct"), "generator"
handlePruned, labelPruned  = Handle ("std::vector<reco::GenParticle>"     ),"prunedGenParticles"
handlePacked, labelPacked  = Handle ("std::vector<pat::PackedGenParticle>"),"packedGenParticles"
fatjets, fatjetLabel       = Handle("std::vector<pat::Jet>"               ),"slimmedJetsAK8"


if infile.fine("/pnfs/")!=-1:
    fname = 'dcap://t3se01.psi.ch:22125/' + infile
else:
    fname = 'root://cms-xrd-global.cern.ch/' + infile

events = Events(fname)

# events = Events('ZprimeToWW_narrow_M-3000_13TeV-madgraph.root')
print "Opening ", fname



for iev,event in enumerate(events):
        if iev % 1000 == 0:print "Event", iev
        # if iev>10: break

        event.getByLabel( genInfoLabel, genInfo )
        genw = genInfo.product().weight()
        
        event.getByLabel( pvLabel, pvInfo )
        PUs = pvInfo.product()
        

        event.getByLabel (labelPruned, handlePruned)
        pruned = handlePruned.product()
        event.getByLabel (labelPacked, handlePacked)
        packed = handlePacked.product()

        # Find hadronic gen W
        IDs = []
        genWs = []
        daughters = []
        
        allGenPs = []
        for p in pruned :
                isLep = False
                daus = []
                if not p.isHardProcess(): continue
                # if p.pt()<1000 or p.pt()>1400 or abs(p.eta())>1.5: continue
                ptlv = ROOT.TLorentzVector()
                ptlv.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass())
                allGenPs.append(ptlv)
                if (isSignal and abs(p.pdgId()) == 24) or (not isSignal and abs(p.pdgId()) != 24) :  
                        for pa in pruned:
                                if not pa.isHardProcess(): continue
                                mother = pa.mother(0)
                                if mother and isAncestor(p,mother):
                                        if 10<abs(pa.pdgId())<19 and pa.isPromptFinalState():
                                                isLep = True
                                        if isLep: continue
                                        if 0<abs(pa.pdgId())<7:
                                            ld = ROOT.TLorentzVector()
                                            ld.SetPtEtaPhiM(pa.pt(),pa.eta(),pa.phi(),pa.mass())
                                            daus.append(ld)
                        if isLep: continue
                        if isSignal and len(daus)==0: continue  
                        l = ROOT.TLorentzVector()
                        l.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),p.mass())
                        genWs.append(l)
                        IDs.append(p.pdgId())
                        daughters.append(daus)
                        
                    
        if isSignal and (len(genWs)==0 or len(daughters)==0): continue

        event.getByLabel(fatjetLabel, fatjets)
        nJ = 0
        for i,j in enumerate(fatjets.product()):
                if j.pt()<200 or abs(j.eta()) > 2.5: continue
                nJ+= 1
                jetID = 0
                # Gen matching
                jtlv = ROOT.TLorentzVector()
                jtlv.SetPtEtaPhiM(j.pt(),j.eta(),j.phi(),j.mass())
                isW = 0
                drMin = 0.6
                index = 99
                for gW in range(0,len(genWs)):
                        
                        dRs = 0
                        if isSignal:
                            for gQs in daughters[gW]:                       
                                    dR = jtlv.DeltaR(gQs)
                                    if dR < 0.8: 
                                        dRs+=1
                            if dRs <2:
                                continue
                            for i in range(2):
                                dauE  [i]=daughters[gW][i].E()  
                                daupt [i]=daughters[gW][i].Pt() 
                                daueta[i]=daughters[gW][i].Eta()
                                dauphi[i]=daughters[gW][i].Phi()

                            vpdgid[0] = IDs[gW]
                            vE[0]     = genWs[gW].E()
                            vpt[0]    = genWs[gW].Pt()
                            veta[0]   = genWs[gW].Eta()
                            vphi[0]   = genWs[gW].Phi()
                            jetID     = IDs[gW]
                            isW=1
                        else:
                            dr = jtlv.DeltaR(genWs[gW])
                            if dr<drMin: 
                                vE[0]     = genWs[gW].E()
                                vpt[0]    = genWs[gW].Pt()
                                veta[0]   = genWs[gW].Eta()
                                vphi[0]   = genWs[gW].Phi()
                                jetID     = IDs[gW]
                                isW = 0
                            else:
                                isW = 1
                                continue
                       
                if isSignal and not isW: continue
                if not isSignal and isW: continue
                
                #Fill tree
                genweight[0] = genw
                i_evt    [0] = iev
                nPV      [0]= PUs[0].getPU_NumInteractions()        
                i_parton [0]  = len(allGenPs)
                type     [0] = isSignal
                jtau21   [0] = j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2")/j.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1")
                je       [0] = j.energy()
                jpt      [0] = j.userFloat("ak8PFJetsPuppiValueMap:pt")
                jeta     [0] = j.userFloat("ak8PFJetsPuppiValueMap:eta")
                jphi     [0] = j.userFloat("ak8PFJetsPuppiValueMap:phi")
                jmass    [0] = j.mass()
                jsd0     [0] = j.userFloat("ak8PFJetsPuppiValueMap:mass")

                chf      [0] = j.chargedHadronEnergyFraction()
                nhf      [0] = j.neutralHadronEnergyFraction()
                nemf     [0] = j.neutralEmEnergyFraction()
                cemf     [0] = j.chargedEmEnergyFraction()
                chMult   [0] = j.chargedMultiplicity()
                neMult   [0] = j.neutralMultiplicity()

                puppi_softdrop, puppi_softdrop_subjet = ROOT.TLorentzVector() ,ROOT.TLorentzVector()
                wSubjets = j.subjets('SoftDropPuppi')
                for iw,wsub in enumerate( wSubjets ) :
                        puppi_softdrop_subjet.SetPtEtaPhiM(wsub.correctedP4(0).pt(),wsub.correctedP4(0).eta(),wsub.correctedP4(0).phi(),wsub.correctedP4(0).mass())
                        puppi_softdrop+=puppi_softdrop_subjet
                jsd0                         [0] = puppi_softdrop.M()

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
                nConst                [0] = len(constituents)
                # lolamatrix = np.zeros((n_features,n_particles))
                
                Es,pxs,pys,pzs  = [np.zeros(n_particles_lola) for _ in range(4)]
                for i2, cand in enumerate(constituents):
                        if i2>(n_particles-1): break
                        if i2<(n_particles_lola):
                            Es [i2]=cand.energy()
                            pxs[i2]=cand.px()    
                            pys[i2]=cand.py()    
                            pzs[i2]=cand.pz()    

                        E[i2]        = cand.energy();
                        x[i2]        = cand.px()    ;
                        y[i2]        = cand.py()    ;
                        z[i2]        = cand.pz()    ;
                        pt[i2]       = cand.pt()
                        vtx[i2]      = cand.vx()
                        vty[i2]      = cand.vy()
                        vtz[i2]      = cand.vz()
                        nPixHits[i2] = cand.numberOfPixelHits()
                        fromPV[i2]   = cand.fromPV()
                        dxy[i2]      = cand.dxy()
                        dxyErr[i2]   = 0.#cand.dxyError()
                        dz[i2]       = cand.dz()
                        dzErr[i2]    = 0.#cand.dzError()
                        charge[i2]   = cand.charge()
                        
                        # fvec = np.array([cand.energy(),cand.px(),cand.py(),cand.pz()])
                        # lolamatrix[:,i2] = fvec
#                         lolamatrix = lolamatrix/500.
                
                
                args = (Es,pxs,pys,pzs)
                df = np.concatenate(args)
                if addLola:
                    X = image_fun(df)
                    # X = image_fun(lolamatrix)
                    score = model.predict_on_batch(X)
                    sigProb =  score[0,1]
                    # print "Score is:" ,score
                    # print "sigProb is:" ,sigProb
                    jlola[0] = sigProb
                                  #
                    # score2 = model_ptW.predict_on_batch(X)
                    # sigProb2 =  score2[0,1]
                    # jlola_ptW[0] = sigProb2
                    
                    # print "Es  " ,Es
                    #print "pxs " ,pxs
                    #print "pys " ,pys
                    #print "pzs " ,pzs
                
                    # df.append(float(isSignal))
                    # if j.userFloat("ak8PFJetsPuppiValueMap:pt") < 1000: continue
                    # df2 = np.append(df,float(isSignal))
                ftuple = tuple(df)
                writer.writerows([ftuple])
                
                t.Fill()


csvf.close()
f.Write()
f.Close()
