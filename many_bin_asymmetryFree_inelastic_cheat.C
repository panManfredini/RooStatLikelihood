
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TLegend.h"
#include "retrieve_input_from_histo.C"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterOriginal.h"
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include <iostream>     // std::cout, std::ios
#include <fstream>     // std::cout, std::ios
#include <stdio.h>

using namespace RooStats ;
using namespace RooFit ;
using namespace std ;

void test_counting_experiment() {

//////////////////////   MODEL BUILDING    /////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/*
N_s = N_tot_theory(Mass,Xsec) * Acceptance_SR * Eff_AmBe_bin_i * mu
N_b = N_Co_SR_bin_i * Norm_factor

Xesec: considered 10^-40 cm^2
Norm_factor = N_Data_CR / N_Co_CR   --> assuming no difference between Co and Data in CR and SR.
N_tot_theory(Mass,Xsec): for 225 livedays, 45kg and considering Xsec. It is a constant, no uncertainty at the moment.

---Costraint Signal
   nuissance parameter = Acceptance_SR, Eff_AmBe_bin_i
   Gauss(Acceptance_SR_obs | Acceptance_SR, err.)
   Poisson(S0_i | S_tot_SR  * Eff_AmBe_bin_i)
   
---Costraint Bkg
   nuissance parameter = N_Co_SR_bin_i, Norm_factor
   Gauss(Norm_factor_obs |  Norm_factor, err)
   Poisson(B0_i | N_Co_SR_bin_i)

---- WARNING:: convergence problems: mu_hat should always be >> 1, too small values have problem in finding minimum
	because mu is set >0. ---> Try to fix Xsec in order to have mu_hat ~ 10

*/

  RooWorkspace w("w");

  //gROOT->ProcessLine(".L retrieve_input_from_histo_NoSys.C+");
  gROOT->ProcessLine(".L retrieve_input_from_histo.C+");

  retrieve_input_from_histo(w);


// Building the model
   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
   mc.SetParametersOfInterest(*w.var("mu"));

// Setting nuissance parameter
   mc.SetNuisanceParameters(*w.set("nuissance_parameter"));

// need now to set the global observable
   mc.SetGlobalObservables(*w.set("g_observables"));

   mc.SetObservables(*w.set("observables"));

// this is needed for the hypothesis tests
   mc.SetSnapshot(*w.var("mu"));


// make data set with the number of observed events
   RooDataSet data("data","", *w.set("observables"));
   data.add(*w.set("observables"));

// import data set in workspace and save it in a file
   w.import(data);

// import model in the workspace 
   w.import(mc);

   w.writeToFile("CountingModel.root", true);




w.Print();

data.Print();

/* 
cout << w.var("S_i")->getValV() << endl;//<< "   "   <<  w.var("S_i_exp")->getValV() << endl;
///////////////////////////////////////////////////////////////////////
ProfileLikelihoodCalculator pl(data,mc);
  pl.SetConfidenceLevel(0.95);
  LikelihoodInterval* interval = pl.GetInterval();

   // find the iterval on the first Parameter of Interest
  RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();

  double lowerLimit = interval->LowerLimit(*firstPOI);
  double upperLimit = interval->UpperLimit(*firstPOI);


  cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
    lowerLimit << ", "<<
    upperLimit <<"] "<<endl;


  LikelihoodIntervalPlot * plot = new LikelihoodIntervalPlot(interval);
//  plot->SetRange(0,50);  // possible eventually to change ranges
  //plot->SetNPoints(50);  // do not use too many points, it could become very slow for some models
  plot->Draw("");  // use option TF1 if too slow (plot.Draw("tf1")

*/




//////////////////////////  hypo test 
  // get the modelConfig (S+B) out of the file
  // and create the B model from the S+B model
  ModelConfig * sbModel = (ModelConfig*) mc.Clone();
  sbModel->SetName("S+B Model");      
  RooRealVar* poi = (RooRealVar*) sbModel->GetParametersOfInterest()->first();
  poi->setVal(1);  // set POI snapshot in S+B model for expected significance
  sbModel->SetSnapshot(*poi);
  ModelConfig * bModel = (ModelConfig*) mc.Clone();
  bModel->SetName("B Model");      
  RooRealVar* poi2 = (RooRealVar*) bModel->GetParametersOfInterest()->first();
  poi2->setVal(0);
  bModel->SetSnapshot( *poi2  );

//------------------Limit calculation for N_th event expected = 10


	  AsymptoticCalculator  ac(data, *bModel, *sbModel);
	  //ac.SetOneSidedDiscovery(true);  // for one-side discovery test
//	  ac.SetOneSided(true);  // for one-side tests (limits)
	    ac.SetQTilde(true);
	  ac.SetPrintLevel(2);  // to suppress print level 


	// create hypotest inverter 
	  // passing the desired calculator 
	  HypoTestInverter *calc = new HypoTestInverter(ac);    // for asymptotic 
	  //HypoTestInverter calc(fc);  // for frequentist

	  calc->SetConfidenceLevel(0.90);
	  //calc->UseCLs(false);
	  calc->UseCLs(true);
	  int npoints = 500;  // number of points to scan
	  //int npoints = 1000;  // number of points to scan default 1000
	  // min and max (better to choose smaller intervals)
	  double poimin = poi->getMin();
	  double poimax = poi->getMax();
	  //poimin = 0; poimax=10;

	  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
	  calc->SetFixedScan(npoints,poimin,poimax);
 	  calc->SetVerbose(2); 
	  HypoTestInverterResult * r = calc->GetInterval();

	  double upperLimit = r->UpperLimit();

	  std::cout << "The computed Expected upper limit is: " <<  r->GetExpectedUpperLimit(0) << std::endl;

//------------ Getting the interval as function of m --------------//
/*   ifstream in;
   in.open("integral_mass.dat");
   

  vector <double> masses_v;
  vector <double> observed_v;
  vector <double> expected_v;
  vector <double> expected_gaud_v;
  vector <double> expected_S1_up_v;
  vector <double> expected_S1_dw_v;
  vector <double> expected_S2_up_v;
  vector <double> expected_S2_dw_v;

  double mass_itr =0.;
  double Nev_exp_th_itr =0.;
  double xsec_modifier = 10.;
  double N_tot_theory = w.var("N_tot_theory")->getValV();

  while(mass_itr <1000.){
	in >> mass_itr;
	in >> Nev_exp_th_itr;

 	
	xsec_modifier = Nev_exp_th_itr * 225.009 * 34.;  //225.009 livedays and 34 kg and 10^-40 cm2 Xsec.

	masses_v.push_back(mass_itr);
	observed_v.push_back( 1.e-40  * N_tot_theory / xsec_modifier * upperLimit );
	expected_v.push_back( 1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(0) );
	expected_gaud_v.push_back(7e-38 *  1.37590955945e-05 / Nev_exp_th_itr );
	expected_S1_up_v.push_back(1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(1));
	expected_S2_up_v.push_back(1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(2));
	expected_S2_dw_v.push_back(1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(-2));
	expected_S1_dw_v.push_back(1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(-1));

	cout << "Expected median limit for mass " << mass_itr << " GeV  = " << 1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(0) << " cm^2 " << endl;
	
//	observed_v.push_back( w.var("Xsec")->getValV() *  w.var("K_m")->getValV()* upperLimit );
//	expected_v.push_back( w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(0) );
//	expected_S1_up_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(1));
//	expected_S2_up_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(2));
//	expected_S2_dw_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(-2));
//	expected_S1_dw_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(-1));


   }


in.close();

const int n = masses_v.size();
double xe[n];
double mA[n];
double observed[n];
double expected[n];
double expected_gaudenz[n];
double exSigma1_l[n];
double exSigma1_u[n];
double exSigma2_l[n];
double exSigma2_u[n];

for(int k=0; k< n; k++){

	mA[k] = masses_v[k];
	observed[k] = observed_v[k];
	expected[k] = expected_v[k];
	expected_gaudenz[k] = expected_gaud_v[k];
	exSigma1_l[k] =expected_v[k] -  expected_S1_dw_v[k] ;
 	exSigma1_u[k] = expected_S1_up_v[k] - expected_v[k];
	exSigma2_l[k] = expected_v[k] - expected_S2_dw_v[k];
	exSigma2_u[k] = expected_S2_up_v[k] - expected_v[k] ;
}

TGraphErrors *obs_limits = new TGraphErrors(n, mA, observed);
TGraphErrors *Exp_limits = new TGraphErrors(n, mA, expected );
TGraphAsymmErrors *Exp_limitsS1 = new TGraphAsymmErrors(n, mA, expected ,xe, xe, exSigma1_l, exSigma1_u );
TGraphAsymmErrors *Exp_limitsS2 = new TGraphAsymmErrors(n, mA, expected ,xe, xe, exSigma2_l, exSigma2_u);

TGraphErrors *Exp_limits_gaudenz = new TGraphErrors( n, mA, expected_gaudenz);

//double expected_xmass[15] = {8e-36,7e-37, 2e-37, 1e-37, 8e-38, 6e-38, 5.5e-38, 5e-38,  4.3e-38, 5e-38, 6e-38, 7e-38, 9e-38, 1.2e-37, 1.5e-37};
//double m_xmass[15] = { 20, 30., 40., 50., 60., 70., 80., 90.,  100., 200., 300., 400, 500.,700., 1000.};

TGraphErrors *Exp_limits_xmass = new TGraphErrors(16);
   Exp_limits_xmass->SetPoint(0,20,8e-36);
   Exp_limits_xmass->SetPointError(0,0,0);
   Exp_limits_xmass->SetPoint(1,29.8071,7.162923e-37);
   Exp_limits_xmass->SetPointError(1,0,0);
   Exp_limits_xmass->SetPoint(2,39.90202,2.027528e-37);
   Exp_limits_xmass->SetPointError(2,0,0);
   Exp_limits_xmass->SetPoint(3,53.41583,9.91722e-38);
   Exp_limits_xmass->SetPointError(3,0,0);
   Exp_limits_xmass->SetPoint(4,62.16429,7.461589e-38);
   Exp_limits_xmass->SetPointError(4,0,0);
   Exp_limits_xmass->SetPoint(5,69.85718,6.3506e-38);
   Exp_limits_xmass->SetPointError(5,0,0);
   Exp_limits_xmass->SetPoint(6,83.21777,5.354015e-38);
   Exp_limits_xmass->SetPointError(6,0,0);
   Exp_limits_xmass->SetPoint(7,90,5e-38);
   Exp_limits_xmass->SetPointError(7,0,0);
   Exp_limits_xmass->SetPoint(8,105.0887,4.600252e-38);
   Exp_limits_xmass->SetPointError(8,0,0);
   Exp_limits_xmass->SetPoint(9,200,5e-38);
   Exp_limits_xmass->SetPointError(9,0,0);
   Exp_limits_xmass->SetPoint(10,300,6e-38);
   Exp_limits_xmass->SetPointError(10,0,0);
   Exp_limits_xmass->SetPoint(11,388.2045,7.252295e-38);
   Exp_limits_xmass->SetPointError(11,0,0);
   Exp_limits_xmass->SetPoint(12,590.8438,9.823615e-38);
   Exp_limits_xmass->SetPointError(12,0,0);
   Exp_limits_xmass->SetPoint(13,746.1269,1.210266e-37);
   Exp_limits_xmass->SetPointError(13,0,0);
   Exp_limits_xmass->SetPoint(14,1000,1.5e-37);
   Exp_limits_xmass->SetPointError(14,0,0);
   Exp_limits_xmass->SetPoint(15,4244.204,4.354065e-37);
   Exp_limits_xmass->SetPointError(15,0,0);


TCanvas *c1 = new TCanvas("limits", "limit", 600, 600);

Exp_limitsS1->SetFillColor(3);
Exp_limitsS1->SetLineColor(3);
Exp_limitsS1->SetMarkerColor(3);
Exp_limitsS1->SetMarkerSize(0);

Exp_limitsS2->SetFillColor(5);
Exp_limitsS2->SetLineColor(5);
Exp_limitsS2->SetMarkerColor(5);
Exp_limitsS2->SetMarkerSize(0);

obs_limits->SetFillColor(0);
obs_limits->SetLineWidth(3);
obs_limits->SetMarkerSize(0);

Exp_limits->SetFillColor(0);
Exp_limits->SetMarkerSize(0);
Exp_limits->SetLineStyle(7);
Exp_limits->SetLineWidth(3);

Exp_limits_gaudenz->SetFillColor(0);
Exp_limits_gaudenz->SetMarkerSize(0);
Exp_limits_gaudenz->SetLineWidth(3);
Exp_limits_gaudenz->SetLineColor(4);

Exp_limits_xmass->SetFillColor(0);
Exp_limits_xmass->SetMarkerSize(0);
Exp_limits_xmass->SetLineWidth(3);
Exp_limits_xmass->SetLineColor(2);

//Exp_limitsS2->GetYaxis()->SetTitle("#sigma#timesBR( #phi #rightarrow #tau#tau )  [pb]");
Exp_limitsS2->GetYaxis()->SetTitle("#sigma");

Exp_limitsS2->GetXaxis()->SetTitle("M  [GeV]");


Exp_limitsS2->GetXaxis()->SetLimits(9.,1000.);
Exp_limitsS2->GetYaxis()->SetRangeUser(1E-38,1E-30);

Exp_limits->GetXaxis()->SetLimits(9.,1000.);
Exp_limits->GetYaxis()->SetRangeUser(1E-38,1E-30);


Exp_limitsS2->Draw("Al3");
Exp_limitsS1->Draw("sameL3");
Exp_limits->Draw("PL");
Exp_limits_gaudenz->Draw("PC");
Exp_limits_xmass->Draw("PC");
//obs_limits->Draw("PL");


TLegend* lego = new TLegend(0.2,0.9,0.5,0.7);
  lego->SetTextSize(0.033);
  lego->SetFillColor(0);
  lego->SetBorderSize(0);
  lego->AddEntry(obs_limits,"Observed 90\% CLs limit");
  lego->AddEntry(Exp_limits_gaudenz, "Expected 90\% Gaudenz");
  lego->AddEntry(Exp_limits_xmass, "Expected 90\% XMASS");
  lego->AddEntry(Exp_limits, "Expected 90\% CLs limit");
  lego->AddEntry(Exp_limitsS1,"1 #sigma","f");
  lego->AddEntry(Exp_limitsS2,"2 #sigma","f");
  lego->Draw();


gPad->SetLogy();
gPad->SetLogx();
gPad->RedrawAxis("g");

myText(0.4,0.86,2,"Test");

*/


  // now use the profile inspector
  ProfileInspector p;
  TList* list = p.GetListOfProfilePlots(data,&mc);
  
  // now make plots
  TCanvas* c1 = new TCanvas("c1","ProfileInspectorDemo",800,200);
  if(list->GetSize()>4){
    double n = list->GetSize();
    int nx = (int)sqrt(n) ;
    int ny = TMath::CeilNint(n/nx);
    nx = TMath::CeilNint( sqrt(n) );
    c1->Divide(ny,nx);
  } else
    c1->Divide(list->GetSize());
  for(int i=0; i<list->GetSize(); ++i){
    c1->cd(i+1);
    list->At(i)->Draw("al");
  }
  
  cout << endl;



/*  // plot now the result of the scan 
  HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot","HypoTest Scan Result",r);

  // plot in a new canvas with style
  TCanvas * c1 = new TCanvas("HypoTestInverter Scan"); 
  c1->SetLogy(false);

  plot->Draw("2CL");  // plot also CLb and CLs+b 
  //plot->Draw("OBS");  // plot only observed p-value

*/



  // plot also in a new canvas the test statistics distributions 
  
  // plot test statistics distributions for the two hypothesis
/*  // when distribution is generated (case of FrequentistCalculators)
  const int n = r->ArraySize();
  if (n> 0 &&  r->GetResult(0)->GetNullDistribution() ) { 
     TCanvas * c2 = new TCanvas("Test Statistic Distributions","",2);
     if (n > 1) {
        int ny = TMath::CeilNint( sqrt(n) );
        int nx = TMath::CeilNint(double(n)/ny);
        c2->Divide( nx,ny);
     }
     for (int i=0; i<n; i++) {
        if (n > 1) c2->cd(i+1);
        SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
        pl->SetLogYaxis(true);
        pl->Draw();
     }
  }
*/




}







  // create the AsymptoticCalculator from data,alt model, null model
/*
 // configure ToyMC Samler (needed only for frequentit calculator)
  ToyMCSampler *toymcs = (ToyMCSampler*)calc.GetHypoTestCalculator()->GetTestStatSampler();
   
  // profile likelihood test statistics 
  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  // for CLs (bounded intervals) use one-sided profile likelihood
  if (useCLs) profll.SetOneSided(true);


   // ratio of profile likelihood - need to pass snapshot for the alt 
  // RatioOfProfiledLikelihoodsTestStat ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   
  // set the test statistic to use 
  toymcs->SetTestStatistic(&profll);
*/
  // if the pdf is not extended (e.g. in the Poisson model) 
  // we need to set the number of events
 // if (!sbModel->GetPdf()->canBeExtended())
   //  toymcs->SetNEventsPerToy(1);





