#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>

// Root
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TMarker.h"
// RooFit
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooHist.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "Roo1DTable.h"
#include "RooConstVar.h"
#include "RooProduct.h"
#include "RooRandom.h"
#include "TStopwatch.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"
#include "RooAbsPdf.h"

// RooStat
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileInspector.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/MinNLLTestStat.h"

using namespace RooStats ;
using namespace RooFit ;
using namespace std ;

void statTest(double mu_pe, double mu_hyp, ModelConfig *mc , RooDataSet *data ){

    int nToyMC = 5;
    // set roofit seed
    RooRandom::randomGenerator()->SetSeed();

    cout << endl;
    cout << endl;
    cout << "Will generate " << nToyMC << " pseudo-experiments for : " << endl;
    cout << " - mu[pseudo-data] = " << mu_pe  << endl;
    cout << " - mu[stat-test]   = " << mu_hyp << endl;
    cout << endl;

    // Check number of POI (for Wald approx)
    RooArgSet *ParamOfInterest = (RooArgSet*) mc->GetParametersOfInterest();
    int nPOI = ParamOfInterest->getSize();
    if(nPOI>1){
      cout <<"not sure what to do with other parameters of interest, but here are their values"<<endl;
      mc->GetParametersOfInterest()->Print("v");
    }
    RooRealVar* firstPOI    = (RooRealVar*) ParamOfInterest->first(); 
    RooAbsPdf *simPdf = (mc->GetPdf());
    //PrintAllParametersAndValues( *mc->GetGlobalObservables() );
    //PrintAllParametersAndValues( *mc->GetObservables() );
    firstPOI->setVal(0.0); // FIXME

    //simPdf->fitTo( *data, Hesse(kTRUE), Minos(kTRUE), PrintLevel(1) );
    simPdf->fitTo( *data );

    // set up the sampler
    ToyMCSampler sampler;
    sampler.SetPdf(*mc->GetPdf());
    sampler.SetObservables(*mc->GetObservables());
    sampler.SetNToys(nToyMC);
    sampler.SetGlobalObservables(*mc->GetGlobalObservables());
    sampler.SetParametersForTestStat(*mc->GetParametersOfInterest());
    RooArgSet* poiset = dynamic_cast<RooArgSet*>(ParamOfInterest->Clone());


    // only unconditional fit
    MinNLLTestStat *minNll = new MinNLLTestStat(*mc->GetPdf());
    minNll->EnableDetailedOutput(true);
    sampler.AddTestStatistic(minNll);

    // enable PROOF if desired
    //ProofConfig pc(*w, 8, "workers=8", kFALSE);
    //sampler.SetProofConfig(&pc);

    // evaluate the test statistics - this is where most of our time will be spent
    cout << "Generating " << nToyMC << " toys...this will take a few minutes" << endl;
    TStopwatch *mn_t = new TStopwatch; 
    mn_t->Start();
    RooDataSet* sd = sampler.GetSamplingDistributions(*poiset);
    cout << "Toy generation complete :" << endl;
    // stop timing
    mn_t->Stop();
    cout << " total CPU time: " << mn_t->CpuTime() << endl;
    cout << " total real time: " << mn_t->RealTime() << endl; 

    // now sd contains all information about our test statistics, including detailed output
    // we might eg. want to explore the results either directly, or first converting to a TTree
    // do the conversion
    TFile f("mytoys.root", "RECREATE");
    TTree *toyTree = RooStats::GetAsTTree("toyTree", "TTree created from test statistics", *sd);
    // save result to file, but in general do whatever you like
    f.cd();
    toyTree->Write();
    f.Close();
/*
    TFile* tmpFile = new TFile("mytoys.root","READ");
    TTree* myTree = (TTree*)tmpFile->Get("toyTree");

    // get boundaries for histograms
    TIter nextLeaf( (myTree->GetListOfLeaves())->MakeIterator() );
    TObject* leafObj(0);
    map<TString, float> xMaxs;
    map<TString, float> xMins;
    for(int i(0); i<myTree->GetEntries(); i++) {
      myTree->GetEntry(i);
      nextLeaf = ( (myTree->GetListOfLeaves())->MakeIterator() );
      while( (leafObj = nextLeaf.Next()) ) {
        TString name(leafObj->GetName());
        float value(myTree->GetLeaf( leafObj->GetName() )->GetValue());
        if(value > xMaxs[name]) { xMaxs[name] = value; }
        if(value < xMins[name]) { xMins[name] = value; }
      } // loop over leaves
    } // loop over tree entries

    // plot everything in the tree
    myTree->GetEntry(0);
    nextLeaf = ( (myTree->GetListOfLeaves())->MakeIterator() );
    leafObj = 0;
    // make a histogram per leaf
    map<TString, TH1F*> hists;
    myTree->GetEntry(0);
    while( (leafObj = nextLeaf.Next()) ) {
      if(!leafObj) { continue; }
      //cout << leafObj->GetName() << endl;
      TString name(leafObj->GetName());
      // special ones : fit related things
      if(name.Contains("covQual"))   { hists[name] = new TH1F(name,name,5,0,5); continue; }
      if(name.Contains("fitStatus")) { hists[name] = new TH1F(name,name,5,0,5); continue; }
      int nbin(500); 
      float histMin( xMins[name] - 0.1*fabs(xMins[name]) ); 
      float histMax( xMaxs[name] + 0.1*fabs(xMaxs[name]) );
      if(name.Contains("ATLAS_norm")) { // floating normalization factors
        histMin = 0; histMax = 10;
      }
      else if(name.Contains("gamma_stat")) { // statistical nus param
        if(name.Contains("globObs")) {  // get custom range for sampling
          histMin = int( xMins[name] - 0.1*fabs(xMins[name]) );
          histMax = int( xMaxs[name] + 0.1*fabs(xMaxs[name]) );
        } // use small range for pull and error
        else { nbin = 100; histMin = 0.0; histMax = 2.0; }
      }
      else if(name.Contains("_err")) { // errors on nus param
        nbin = 100; histMin = 0.0; histMax = 2.0;
      }
      else if(name.Contains("fitCond") || name.Contains("fitUncond") || name.Contains("globObs")) { // fit pulls
        nbin = 500; histMin = -5; histMax = 5;
      }
      hists[name] = new TH1F(name,name,nbin,histMin,histMax);
    } // loop over leaves to declare histos

    // loop over entries and fill histograms
    for(int i(0); i<myTree->GetEntries(); i++) {
      myTree->GetEntry(i);
      nextLeaf = ( (myTree->GetListOfLeaves())->MakeIterator() );
      while( (leafObj = nextLeaf.Next()) ) {
        TString name(leafObj->GetName());
        if(hists.find(name) == hists.end()) { continue; }
        hists[name]->Fill( myTree->GetLeaf( leafObj->GetName() )->GetValue() );
      } // loop over leaves
    } // loop over tree entries

    // overflow and underflow
    for(map<TString,TH1F*>::iterator ihist(hists.begin()); ihist!=hists.end(); ihist++) {
      if(ihist->second->GetBinContent(0)>0) {
        ihist->second->SetBinContent(1, ihist->second->GetBinContent(0) + ihist->second->GetBinContent(1) );
        // fix err
      }
      int nBinx = ihist->second->GetNbinsX();
      if(ihist->second->GetBinContent(nBinx)>0) {
        ihist->second->SetBinContent(nBinx-1, ihist->second->GetBinContent(nBinx) + ihist->second->GetBinContent(nBinx-1) );
        // fix err
      }
    }

    // save the results
    TString dirName(OutputDir+"/PlotsStatisticalTest/GlobalFit");
    if(drawPlots) {
      system(TString("mkdir -vp "+dirName));
    }
    TCanvas* canvas = new TCanvas("pulls");
    TLegend *leg = new TLegend(0.67, 0.64, 0.87, 0.86);
    LegendStyle(leg);
    for(map<TString,TH1F*>::iterator ihist(hists.begin()); ihist!=hists.end(); ihist++) {
      if( (ihist->first).Contains("fitCond_") ) { continue; } // skip unconditional fit - get it explicitly
      canvas->Clear();
      leg->Clear();
      TString niceName(ihist->first);
      niceName.ReplaceAll("fitUncond_","");
      //niceName.ReplaceAll("SD_TS0_",""); // not good if have multiple test statistics
      // conditional fit information
      ihist->second->SetLineColor(kGray+2);
      ihist->second->SetTitle(niceName);
      ihist->second->SetLineStyle(kSolid);
      ihist->second->SetLineWidth(2);
      if((ihist->first).Contains("fit") && !(ihist->first).Contains("_err") 
          && !(ihist->first).Contains("Qual") && !(ihist->first).Contains("Status")) {
        ihist->second->Rebin(4);
      }

//      ihist->second->GetXaxis()->SetTitle("");
//      ihist->second->GetYaxis()->SetTitle("");

      if(niceName.Contains("globObs")) {
        leg->AddEntry( ihist->second, "Sampling", "l" ); // add value of mu
      } else {
        leg->AddEntry( ihist->second, "Unconditional Fit", "l" ); // add value of mu
      }
      TString condName(ihist->first);
      condName.ReplaceAll("fitUncond","fitCond");
      // uncomditional fit information
      if(hists.find(condName) != hists.end() && condName != ihist->first) {
        hists[condName]->SetLineColor(kGray+2);
        hists[condName]->SetLineStyle(kDashed);
        hists[condName]->SetLineWidth(2);
        if(!(ihist->first).Contains("_err")) { hists[condName]->Rebin(4); }
        leg->AddEntry( hists[condName], "Conditional Fit", "l" );
        if( hists[condName]->GetMaximum() > ihist->second->GetMaximum() ) {
          ihist->second->SetMaximum( hists[condName]->GetMaximum() );
        }
      }
      ihist->second->SetMaximum( 1.2 * ihist->second->GetMaximum() );
      canvas->cd();
      ihist->second->Draw();
      leg->Draw();
      if(hists[condName] && condName != ihist->first) { hists[condName]->Draw("same"); }
      if(drawPlots) { 
        canvas->Print(dirName+"/"+niceName+".eps");
        canvas->Print(dirName+"/"+niceName+".png");
      }

      MainDirStatTest->cd();
      canvas->Write();
      gROOT->cd();
    }

*/

return;
}
