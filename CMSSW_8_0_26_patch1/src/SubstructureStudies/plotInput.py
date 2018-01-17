from ROOT import TFile, TH1F, TCanvas, TTree, TLegend, TGraph, TColor
from time import sleep
import tdrstyle
from CMS_lumi import *
tdrstyle.setTDRStyle()
import ROOT 

ROOT.gROOT.SetBatch(True)

H_ref = 600
W_ref = 800
W = W_ref
H  = H_ref
T = 0.08*H_ref
B = 0.12*H_ref 
L = 0.12*W_ref
R = 0.04*W_ref

CMS_lumi.lumi_13TeV = ""
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 0
lumi = "1"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
	
def getPavetext():
  addInfo = TPaveText(0.3010112,0.2066292,0.4202143,0.3523546,"NDC")
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  return addInfo
  
def get_palette(mode):
 palette = {}
 palette['gv'] = [] 
 colors = ['#FF420E','#80BD9E','#336B87','#763626','#003B46','#66A5AD','#2E4600']
 for c in colors:
	 palette['gv'].append(c)
 return palette[mode]
 
palette = get_palette('gv')
col = TColor() 

def drawTH1(tree,var,cuts,bins,min,max,fillcolor,titlex = "",units = "",drawStyle = "HIST"):
	h = ROOT.TH1D("tmpTH1","",bins,min,max)
	h.Sumw2()
	# h.SetLineColor(1)
# 	h.SetLineWidth(2)
# 	h.SetFillStyle(1)
	h.SetFillColor(fillcolor)
	h.GetYaxis().SetTitle("A.U")
	h.GetYaxis().SetTitleOffset(1.0)
	h.GetYaxis().SetNdivisions(414)
	if units=="":
	    h.GetXaxis().SetTitle(titlex)
	else:
	    h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
	corrString='1'
	tree.Draw(var+">>tmpTH1","("+cuts+")*"+lumi+"*("+corrString+")","goff")
	return h
	


vars=["ve","vpt","vpx","vpy","vpz","nconst","jflavour","vmass","msoftdrop_beta0","msoftdrop_beta2","je","jpt","jeta","jphi","jtau21","jtau21sd0","jtau21sd2","pe","ppt","ppx","ppy","ppz"]
titlex=["W energy","W p_{t}","W p_{x}","W p_{y}","W p_{z}","Number of constrituents","Jet flavour","W mass","AK8 mMDT","AK8 SD #beta=2","AK8 energy","AK8 p_{T}","AK8 #eta","AK8 #phi","AK8 #tau_{21}","AK8 #tau_{21}^{mMDT}","AK8 #tau_{21}^{SD #beta=2}","Jet constituent energy","Jet constituent p_{T}","Jet constituent p_{x}","Jet constituent p_{y}","Jet constituent p_{z}"]
unit=["GeV","GeV","GeV","GeV","GeV","","","GeV","GeV","GeV","GeV","GeV","","","","","","GeV","GeV","GeV","GeV","GeV"]
minx=[0,0,-200,-200,-200,0,0,0,0,0,0,0,-3.2,-3.2,0,0,0,0,0,-10,-10,-10]
maxx=[3700,3700,200,200,200,200,3,300,300,300,3700,3700,3.2,3.2,1,1,1,10,10,10,10,10]
binsx=[74,74,40,40,40,200,3,60,60,60,74,74,40,40,20,20,20,40,40,40,40,40]

# vars=["ve"]

files = []
sf = TFile.Open("/scratch/thaarres/SUBSTRUCTURE/LOLAoutput/Signal.root","READ")
bf = TFile.Open("/scratch/thaarres/SUBSTRUCTURE/LOLAoutput/Background.root","READ")
files.append(sf)
files.append(bf)
for i,var in enumerate(vars):
	

	name = var
	canvas = ROOT.TCanvas(name+"_c1",name+"_c1",50,50,W,H)
	canvas.SetFillColor(0)
	canvas.SetBorderMode(0)
	canvas.SetFrameFillStyle(0)
	canvas.SetFrameBorderMode(0)
	canvas.SetLeftMargin( L/W )
	canvas.SetRightMargin( R/W )
	canvas.SetTopMargin( T/H )
	canvas.SetBottomMargin( B/H )
	canvas.SetTickx(0)
	canvas.SetTicky(0)
	
	
	
	legend = ROOT.TLegend(0.62,0.7,0.92,0.9,"","brNDC")
	legend.SetBorderSize(0)
	legend.SetLineColor(1)
	legend.SetLineStyle(1)
	legend.SetLineWidth(1)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	
	hists=[]
	style = [1,3002]
	for ii,file in enumerate(files):	
		tree = file.Get('tree')
		cutL="1"
		hist = drawTH1(tree,var,cutL,binsx[i],minx[i],maxx[i],col.GetColor(palette[ii]),titlex[i],unit[i])
		hist.SetName(name+str(i))
		ROOT.SetOwnership(hist,True)
		print hist.GetName()
		if file.GetName().find("Signal")!=-1: legend.AddEntry(hist,"W (G_{Bulk}#rightarrowWW)","F")
		else: legend.AddEntry(hist,"QCD (Herwig++)","F")
		hist.SetFillColor(col.GetColor(palette[ii]))
		hist.SetFillStyle(style[ii])
		hist.Scale(1./hist.Integral())
		hists.append(hist)
	

	canvas.cd()
	
	hists[0].Draw("HIST")
	maxy = hists[0].GetMaximum()
	if maxy < hists[1].GetMaximum():
		maxy = hists[1].GetMaximum()
	hists[0].GetYaxis().SetRangeUser(0, maxy*1.5);
	for h in hists:
		h.Draw("sameHIST")


	legend.Draw("SAME")
	CMS_lumi(canvas, iPeriod, iPos)
	canvas.Update()
	# sleep(100)
	canvas.SaveAs("controlplots/"+var+".png")
sf.Close()
bf.Close()
del files
del hists
del canvas
	

	