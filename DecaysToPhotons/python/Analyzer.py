from DataFormats.FWLite import Events, Handle
from CMSOpenDataAnalysis.DecaysToPhotons import rootnotes
import ROOT
import os
import json
import time

def getFullPath(path):
    return os.path.join(os.environ['CMSSW_BASE'],
                        'src/CMSOpenDataAnalysis/DecaysToPhotons',
                        path)


class EventBox(object):
    def __init__(self):
        self.pho1 = None
        self.pho2 = None


class Object(object):
    """
    What is this?
    """
    def __init__(self, pho1, pho2=None):
        self.pho1 = pho1
        self.pho2 = pho2
        if pho2:
            self.p4 = ROOT.TLorentzVector(self.pho1.px() + self.pho2.px(),
                                          self.pho1.py() + self.pho2.py(),
                                          self.pho1.pz() + self.pho2.pz(),
                                          self.pho1.energy() + self.pho2.energy())
        else:
            self.p4 = ROOT.TLorentzVector(self.pho1.px(),
                                          self.pho1.py(),
                                          self.pho1.pz(),
                                          self.pho1.energy())

    def mass(self):
        return self.p4.M()

#    def pdgId(self):
#        return 25

    def p4(self):
        return self.p4

    def pt(self):
        return self.p4.Pt()

    def px(self):
        return self.p4.Px()

    def py(self):
        return self.p4.Py()

    def pz(self):
        return self.p4.Pz()

    def energy(self):
        return self.p4.Energy()

    def eta(self):
        return self.p4.eta()

    def phi(self):
        return self.p4.Phi()

    def rapidity(self):
        return self.p4.Rapidity()

class Analyzer (object):
    """
    Skeleton for building Analyzers
    """

    def __init__(self):
#        self.vertexHandle = Handle('std::vector<reco::Vertex>')
        self.photonHandle = Handle('std::vector<reco::Photon>')
        self.histograms = {}
        self.samples = ['data']
        self.plots = {}
        for sample in self.samples:
            self.histograms[sample] = {}
        self.data = []

#    def photonID(self, photon, vertex):
    def photonID(self, photon):
        """
        Subclasses of Analyzer should define this method
        """
        return True


    def readCollections(self, event, box, isFake=False):
#        event.getByLabel('offlinePrimaryVertices', self.vertexHandle)
        event.getByLabel('photons', self.photonHandle)

        box.photons = self.photonHandle.product()
#        box.vertex = self.vertexHandle.product()[0]

        #first select photons
        box.selectedPhotons = []
        #PT>10. and |eta|<2.5 but not the EB/EE gap 
        for pho in box.photons:
 	    #if pho.pt() < 10 or abs(pho.eta()) > 2.5 or (abs(pho.eta()) > 1.4442 and abs(pho.eta()) < 1.566):
	    SCeta=pho.superCluster().position().Eta()
            if pho.pt() < 10 or abs(SCeta) > 2.5 or (abs(SCeta) > 1.4442 and abs(SCeta) < 1.566):
                continue
            try:
#                if self.photonID(pho, box.vertex):
                if self.photonID(pho):
                    box.selectedPhotons.append(pho)
            except ZeroDivisionError:
                continue

    def exportData(self):
        f = open('data.json', 'w')
        info = {'data': self.data}
        json.dump(info, f)
        f.close()

    def analyze(self, box):
        return True

    def declareHistos(self):
        return True

    def fillHistos(self, box, sample, weight=1):
        return True

    def declareHisto(self, name, bins, min, max, xlabel=''):
        for sample in self.samples:
            self.histograms[sample][name] = ROOT.TH1F(sample + '_' + name,
                                                      name, bins, min, max)
            self.histograms[sample][name].GetXaxis().SetTitle(xlabel)
            self.histograms[sample][name].GetYaxis().SetTitle('events')

    def fillHisto(self, name, sample, value, weight=1):
        self.histograms[sample][name].Fill(value, weight)

    def addEvent(self, box):
        return

    def processSample(self, sample, maxEv=-1):
        print 'Processing Files'
        print sample.files

        events = Events(sample.files)
        print "%s events available for processing" % events.size()
        ts = time.time() 
        for N, event in enumerate(events):
            if maxEv >= 0 and (N + 1) >= maxEv:
               break
            if N % 1000000 == 0:
                t2 = time.time()
                print "%s events processed in %s seconds" % (N + 1, t2 - ts)
            weight = 1
            box = EventBox()
            self.readCollections(event, box)
            if not self.analyze(box):
                continue
            self.addEvent(box)
            self.fillHistos(box, sample.type, weight)
        tf = time.time()
        print "%s events processed in %s seconds" % (N + 1, tf - ts)

    def convertToPoisson(self, h):
        graph = ROOT.TGraphAsymmErrors()
        q = (1 - 0.6827) / 2.

        for i in range(1, h.GetNbinsX() + 1):
            x = h.GetXaxis().GetBinCenter(i)
            # XXX NOT USED!
            # xLow = h.GetXaxis().GetBinLowEdge(i)
            # xHigh = h.GetXaxis().GetBinUpEdge(i)
            y = h.GetBinContent(i)
            yLow = 0
            yHigh = 0
            if y != 0.0:
                yLow = y - ROOT.Math.chisquared_quantile_c(1 - q, 2 * y) / 2.
                yHigh = ROOT.Math.chisquared_quantile_c(q,
                                                        2 * (y + 1)) / 2. - y
                graph.SetPoint(i - 1, x, y)
                graph.SetPointEYlow(i - 1, yLow)
                graph.SetPointEYhigh(i - 1, yHigh)
                graph.SetPointEXlow(i - 1, 0.0)
                graph.SetPointEXhigh(i - 1, 0.0)

        graph.SetMarkerStyle(20)
        graph.SetLineWidth(2)
        graph.SetMarkerSize(1.)
        graph.SetMarkerColor(ROOT.kBlack)

        return graph

    def makePlot(self, histogram):
        # XXX NOT USED!
        # sandbox = []
        #canvas = ROOT.TCanvas(histogram)
        #if histogram in self.plots:
        #    return self.plots[histogram]['canvas']
        canvas = rootnotes.canvas(histogram)
        canvas.cd()
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        canvas.Range(-68.75, -7.5, 856.25, 42.5)
        canvas.SetFillColor(0)
        canvas.SetBorderMode(0)
        canvas.SetBorderSize(2)
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.SetLeftMargin(0.15)
        canvas.SetRightMargin(0.05)
        canvas.SetTopMargin(0.05)
        canvas.SetBottomMargin(0.15)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.cd()

        frame = canvas.DrawFrame(
            self.histograms['data'][histogram].GetXaxis().GetXmin(),
            0.0, self.histograms['data'][histogram].GetXaxis().GetXmax(), 
            self.histograms['data'][histogram].GetMaximum()*1.2)

        frame.GetXaxis().SetLabelFont(42)
        frame.GetXaxis().SetLabelOffset(0.007)
        frame.GetXaxis().SetLabelSize(0.045)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetTitleOffset(1.15)
        frame.GetXaxis().SetTitleFont(42)
        frame.GetYaxis().SetLabelFont(42)
        frame.GetYaxis().SetLabelOffset(0.007)
        frame.GetYaxis().SetLabelSize(0.045)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleOffset(1.4)
        frame.GetYaxis().SetTitleFont(42)
        frame.GetZaxis().SetLabelFont(42)
        frame.GetZaxis().SetLabelOffset(0.007)
        frame.GetZaxis().SetLabelSize(0.045)
        frame.GetZaxis().SetTitleSize(0.05)
        frame.GetZaxis().SetTitleFont(42)

        frame.GetXaxis().SetTitle(self.histograms['data'][histogram].GetXaxis().GetTitle())
        frame.GetYaxis().SetTitle(self.histograms['data'][histogram].GetYaxis().GetTitle())
        frame.Draw()

        self.histograms['data'][histogram].SetMarkerStyle(20)
        self.histograms['data'][histogram].Sumw2()
        dataG = None
        if self.histograms['data'][histogram].Integral() > 0:
            dataG = self.convertToPoisson(self.histograms['data'][histogram])
            dataG.Draw("Psame")

        legend = ROOT.TLegend(0.74, 0.84, 0.94, 0.94, "", "brNDC")
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.AddEntry(self.histograms['data'][histogram], "Data", "p")
        legend.Draw()

        pt = ROOT.TPaveText(0.1577181, 0.9562937, 0.9580537,
                            0.9947552, "brNDC")
        pt.SetBorderSize(0)
        pt.SetTextAlign(12)
        pt.SetFillStyle(0)
        pt.SetTextFont(42)
        pt.SetTextSize(0.03)
        pt.AddText(0.01, 0.4, "CMS")
        pt.Draw()

        self.plots[histogram] = {
            'canvas': canvas,
            'legend': legend,
            'dataG': dataG,
            'latex1': pt}

        canvas.RedrawAxis()
        canvas.Update()

        canvas.SaveAs(histogram+".png")          
       
        return canvas

    def makeAllPlots(self):
        return [self.makePlot(h) for h in self.histograms['data']]

##############################################################################
    def saveAllHists(self, OutputFileName):

        file=ROOT.TFile(OutputFileName,"RECREATE")
        file.cd()
        for h in self.histograms['data']:
            self.histograms['data'][h].Write()
            
        file.Write()
        file.Close()
        return 
