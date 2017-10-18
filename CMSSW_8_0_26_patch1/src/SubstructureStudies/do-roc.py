from ROOT import TFile, TH1F, TCanvas, TTree, TLegend, TGraph, TColor
from time import sleep
import numpy as np
import h5py
import tdrstyle,CMS_lumi
tdrstyle.setTDRStyle()

CMS_lumi.lumi_13TeV = ""
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Simulation Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

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

def get_stdroc(Signal, Background):
	
	sf = TFile.Open(Signal		+".root","READ")
	bf = TFile.Open(Background	+".root","READ")
	hists = ["h_tau21_mmdt","h_tau21_sd","h_tau21_dichr","h_tau21_puppi"]#,"h_n2_mmdt","h_n2_sd","h_n2_dichr"]
	
	
	
	S_sd = sf.Get("h_sd")
	B_sd = bf.Get("h_sd")
	S_puppi = sf.Get("h_puppimass")
	B_puppi = bf.Get("h_puppimass")
	S_sd.SetLineColor(col.GetColor(palette[0]))
	B_sd.SetLineColor(col.GetColor(palette[0]))
	S_puppi.SetLineColor(col.GetColor(palette[1]))
	B_puppi.SetLineColor(col.GetColor(palette[1]))
	
	histsCut = []
	leg = TLegend(0.65,0.65,0.88,0.85)
	colors = [1,2,3]
	
	#Create histograms
	ShistDEN = sf.Get("h_mmdt")
	ShistDEN.SetName("S_mmdt")
	BhistDEN = bf.Get("h_mmdt")
	ShistDEN.SetName("B_mmdt")

	ShistDEN.SetLineColor(col.GetColor(palette[2]))
	BhistDEN.SetLineColor(col.GetColor(palette[2]))
	ShistsCut = []
	BhistsCut = []
	
	for i,h in enumerate(hists):
		print "Hist number " ,h
		histSCut = sf.Get(h)
		histSCut.SetLineColor(col.GetColor(palette[i]))
		histSCut.SetLineWidth(3)
		histSCut.SetName("S"+h)
		ShistsCut.append(histSCut)
		leg.AddEntry(histSCut,histSCut.GetName().replace("Sh_tau21_"," "),"L")
		
		histBCut = bf.Get(h)
		histBCut.SetLineColor(col.GetColor(palette[i]))
		histBCut.SetLineStyle(2)
		histBCut.SetLineWidth(3)
		histBCut.SetName("B"+h)
		BhistsCut.append(histBCut)
		
	
	
	leg2 = TLegend(0.73,0.3,0.89,0.47)
	leg2.AddEntry(ShistsCut[2],"W'#rightarrowWZ","L")	
	leg2.AddEntry(BhistsCut[2],"QCD","L")
	
	leg = TLegend(0.63,0.17,0.93,0.48)
	name = ["#tau_{21}^{mmdt}+m_{mMDT}","#tau_{21}^{sd}+m_{mMDT}","#tau_{21}^{dichr}+m_{mMDT}","#tau_{21}^{PUPPI}+m_{mMDT}","N_{2}^{mmdt}+m_{mMDT}","N_{2}^{sd}+m_{mMDT}","N_{2}^{dichr}+m_{mMDT}"]
	graphs = []
	
	denom_S = float( ShistDEN.Integral(0, ShistDEN.GetNbinsX()) )
	denom_B = float( BhistDEN.Integral(0, BhistDEN.GetNbinsX()) )
	i = -1
	for hS,hB in zip(ShistsCut,BhistsCut):
		sarray = []
		barray = []
		i += 1
		wpMass  = TGraph(1)
		nbins  	= hS.GetNbinsX()
		g 		= TGraph(nbins)
		
		for ii in range(nbins):
			num_S = float( hS.Integral(0,ii) ) ##Cut from min bin to higher
			num_B = float( hB.Integral(0,ii) )
			g.SetPoint( ii, (num_S/denom_S), (num_B/denom_B) )     #eS vs. eB
			sarray.append(num_S/denom_S)
			barray.append(num_B/denom_B)
		graphs.append(g)
		snparray = np.array(sarray)
		bnparray = np.array(barray)
	
		h5f = h5py.File(hS.GetName()+'std_roc.h5', 'w')
		h5f.create_dataset('dataset_s', data=snparray)
		h5f.create_dataset('dataset_b', data=bnparray)
		h5f.close()
		
		leg.AddEntry(g,name[i],"L")	
	cr = TCanvas()
	cr.SetLogy()
	graphs[0].GetXaxis().SetRangeUser(0.0,1.0)
	graphs[0].GetYaxis().SetRangeUser(0.001,1.0)
	graphs[0].GetYaxis().SetNdivisions(507)
	graphs[0].GetXaxis().SetNdivisions(507)
	graphs[0].Draw("AL")
	
	for i,g in enumerate(graphs):
		g.GetXaxis().SetTitle("Signal efficiency")
		g.GetYaxis().SetTitle("Mistagging rate")
		
		g.SetLineColor(col.GetColor(palette[i]))
		g.SetLineWidth(2)
		g.Draw("Lsame")
	leg.Draw("same")
	CMS_lumi.CMS_lumi(cr, iPeriod, iPos)
	cr.SaveAs("roc.png")
	
	
	
	
	c1 = TCanvas()
	leg = TLegend(0.670,0.60,0.89,0.90)
	i = -1
	for hS,hB in zip(ShistsCut,BhistsCut):
		i += 1
		leg.AddEntry(hS,name[i],"L")	
		hS.Rebin(5)
		hB.Rebin(5)
		hS.GetYaxis().SetRangeUser(0,hS.GetMaximum()*1.5)
		hS.GetXaxis().SetTitle("#tau_{21}")
		hB.GetXaxis().SetTitle("#tau_{21}")
		hS.GetYaxis().SetTitle("A.U")
		hB.GetYaxis().SetTitle("A.U")
		hS.DrawNormalized("HISTCsame")
		hS.GetYaxis().SetNdivisions(507)
		hS.GetXaxis().SetNdivisions(507)
		hB.DrawNormalized("HISTCsame")
	
	leg.Draw("same")	
	leg2.Draw("same")
	CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
	c1.SaveAs("tau21plots.png")
	
	leg = TLegend(0.65,0.65,0.90,0.85)
	leg.AddEntry(ShistDEN,"PUPPI mMDT","L")	
	leg.AddEntry(S_sd    ,"PUPPI SD_{#beta=2}","L")
	leg.AddEntry(S_puppi ,"PUPPI plain mass","L")
	
	plotHists = []
	plotHists.append(ShistDEN)
	plotHists.append(BhistDEN)
	plotHists.append(S_sd	)
	plotHists.append(B_sd	)
	plotHists.append(S_puppi	)
	plotHists.append(B_puppi	)
	
	
	c = TCanvas()
	c.cd()
	for j,h in enumerate(plotHists): 
		h.Rebin(5)
		
		if j%2 ==0: h.SetLineStyle(1)
		else: h.SetLineStyle(2)
		h.SetLineWidth(3)
		h.GetXaxis().SetTitle("Mass (GeV)")
		h.GetYaxis().SetTitle("A.U")
		
		h.DrawNormalized("CHISTsame")
		h.GetYaxis().SetNdivisions(507)
		h.GetXaxis().SetNdivisions(507)
		
	leg.Draw("same")	
	leg2.Draw("same")	
	CMS_lumi.CMS_lumi(c, iPeriod, iPos)
	c.SaveAs("massplots.png")
	
	
	sleep(100)
	return snparray,bnparray
		
if __name__ == "__main__":
	snarray,bnarray = get_stdroc("dichroic_sig","dichroic_bkg")
	