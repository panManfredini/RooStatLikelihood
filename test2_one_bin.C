
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TLegend.h"

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

double N_bkg_SR = 24318.;
double n_bins = 1.;

TString total_model = "PROD:model(";


w.factory("Acc[0.9]");
w.factory("eff[0.5]");
w.factory("mu[1.,0.,100]");
w.factory("Bkg_norm_factor[0.03,0,1]");
w.factory("N_tot_theory[1]");

w.var("eff")->setVal(1./n_bins);
w.var("Acc")->setConstant(true); 
w.var("eff")->setConstant(true); 

RooArgSet observables;
RooArgSet g_obs;
RooArgSet nuissanceP;

w.factory("Gaussian:costraint_bkg_norm(Bkg_norm_factor_obs[0,1], Bkg_norm_factor, err_Norm_factor[0.01])"); // approximation (should be kind of complicated cauchy....)

w.var("N_tot_theory")->setVal(1); 		// Total expected Signal events
w.var("Bkg_norm_factor_obs")->setVal(0.0378); // Ratio between data in CR and Co in CR
w.var("err_Norm_factor")->setVal(0.00013); 		/// Set error on norm factor

w.var("N_tot_theory")->setConstant(true); 
w.var("Bkg_norm_factor_obs")->setConstant(true); 
w.var("err_Norm_factor")->setConstant(true); 

w.var("Bkg_norm_factor")->setMax(w.var("Bkg_norm_factor_obs")->getValV() + w.var("err_Norm_factor")->getValV()*3.); //set to 10 sigma.
w.var("Bkg_norm_factor")->setVal(w.var("Bkg_norm_factor_obs")->getValV());// Set to observed value!!!


   w.factory("prod:N_s(N_tot_theory, Acc, eff, mu)");
   w.factory("prod:N_b(N_Co_SR_bin[1000,0,50000], Bkg_norm_factor)");
   w.factory("sum:nexp(N_s, N_b)");

   w.factory("Poisson:pdf(nobs[0,1000],nexp)");

   //constraint, stat uncertainty
//   w.factory("Gamma::gamma(N_Co_SR_bin,sum::temp(Cobs[1000],1),1,0)") ;
//   w.factory("Poisson:Co_constrain(Cobs[50000],N_Co_SR_bin)");
   
   w.factory("err[5]");
   w.var("err")->setVal(0.0005);
   w.var("err")->setConstant(true);
	
   w.factory("Gaussian:Co_constrain(Cobs[50000],N_Co_SR_bin,err)");


   w.var("nobs")->setVal(N_bkg_SR*0.0378/n_bins);
   w.var("N_Co_SR_bin")->setVal(N_bkg_SR/n_bins);
   w.var("Cobs")->setVal(N_bkg_SR/n_bins);
   w.var("Cobs")->setConstant(true); 

//  w.var("N_Co_SR_bin")->setConstant(true);  

  // total_model.Append("pdf,gamma,");

   observables.add(*w.var("nobs"));

   g_obs.add(*w.var("Cobs"));
   nuissanceP.add(*w.var("N_Co_SR_bin"));


//------------priors   
//  w.factory("Uniform::prior_nuis({Bkg_norm_factor,N_Co_SR_bin_1})");
//  w.factory("Uniform::prior_nuis({Bkg_norm_factor,N_Co_SR_bin_1,N_Co_SR_bin_2})");
//  w.factory("Uniform::prior_poi({mu})");
//  w.factory("PROD::prior(prior_poi,prior_nuis)"); 

w.factory("PROD::model(pdf,Co_constrain,costraint_bkg_norm)");

w.Print();

/*// Variable definition
   w.factory("prod:N_s_1(N_tot_theory[5], Acc[0.9], eff[0.5], mu[1.,0.,100])");
   w.factory("prod:N_b_1(N_Co_SR_bin_1[12159], Bkg_norm_factor[0.03, 0.0,1.0])"); 
   w.factory("sum:nexp_1(N_s_1, N_b_1)");


   w.factory("prod:N_s_2(N_tot_theory[5], Acc[0.9], eff[0.5], mu[1.,0.,100])");
   w.factory("prod:N_b_2(N_Co_SR_bin_2[12159], Bkg_norm_factor[0.03, 0.0,1.0])"); 
   w.factory("sum:nexp_2(N_s_1, N_b_1)");

// Poisson of (n | s+b)
   w.factory("Poisson:pdf_1(nobs_1[100,1000],nexp_1)");
   w.factory("Poisson:pdf_2(nobs_2[100,1000],nexp_2)");


// Constraint of GLOBAL nuissance parameter
   w.factory("Gaussian:costraint_bkg_norm(Bkg_norm_factor_obs[0,1], Bkg_norm_factor, err_Norm_factor[0.01])"); // approximation (should be kind of complicated cauchy....)

// the total model p.d.f will be the product of the two
   w.factory("PROD:model(pdf_1, pdf_2, costraint_bkg_norm)");

// Set input variables -- Global OBSERVABLESerr_Norm_factor
   w.var("N_tot_theory")->setVal(1); 		// Total expected Signal events
   w.var("Bkg_norm_factor_obs")->setVal(0.0378); // Ratio between data in CR and Co in CR
   w.var("err_Norm_factor")->setVal(0.00013); 		/// Set error on norm factor

   

   w.var("N_tot_theory")->setConstant(true); 
   w.var("Acc")->setConstant(true); 
   w.var("eff")->setConstant(true); 
   w.var("N_Co_SR_bin_1")->setConstant(true); 
   w.var("N_Co_SR_bin_2")->setConstant(true); 
   w.var("Bkg_norm_factor_obs")->setConstant(true); 
   w.var("err_Norm_factor")->setConstant(true); 

// Set input Observables
  w.var("nobs_1")->setVal(0.0378*24318./2.);
  w.var("nobs_2")->setVal(0.0378*24318./2.);

// Set range and value for nuissance parameters
   w.var("Bkg_norm_factor")->setMax(w.var("Bkg_norm_factor_obs")->getValV() + w.var("err_Norm_factor")->getValV()*3.); //set to 10 sigma.
   w.var("Bkg_norm_factor")->setVal(w.var("Bkg_norm_factor_obs")->getValV());// Set to observed value!!!

*/


// Building the model
   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
   mc.SetParametersOfInterest(*w.var("mu"));
   //mc.SetPriorPdf(*w.pdf("prior_s"));

// Setting nuissance parameter
   		//RooArgSet nuissanceParameter;  // No Sys
   nuissanceP.add(*w.var("Bkg_norm_factor"));	
   mc.SetNuisanceParameters(nuissanceP);

//   mc.SetPriorPdf(*w.pdf("prior"));
// need now to set the global observable
  		 //RooArgSet g_obs;  // noSYS
   g_obs.add(*w.var("Bkg_norm_factor_obs"));
   mc.SetGlobalObservables(g_obs);

   mc.SetObservables(observables);


// this is needed for the hypothesis tests
   mc.SetSnapshot(*w.var("mu"));


// make data set with the number of observed events
RooDataSet data("data","", observables);
data.add(observables);

// import data set in workspace and save it in a file
   w.import(data);

// import model in the workspace 
   w.import(mc);

 //  w.writeToFile("CountingModel.root", true);




w.Print();

//data.Print();





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

//------------------Limit calculation for N_th event expected = 1


	  AsymptoticCalculator  ac(data, *bModel, *sbModel);
	  //ac.SetOneSidedDiscovery(true);  // for one-side discovery test
	  ac.SetOneSided(true);  // for one-side tests (limits)
	  //  ac->SetQTilde(true);
	  ac.SetPrintLevel(2);  // to suppress print level 

	  

	// create hypotest inverter 
	  // passing the desired calculator 
	  HypoTestInverter *calc = new HypoTestInverter(ac);    // for asymptotic 
	  //HypoTestInverter calc(fc);  // for frequentist


	  calc->SetConfidenceLevel(0.90);
	  calc->UseCLs(true);
	  int npoints = 1000;  // number of points to scan
	  // min and max (better to choose smaller intervals)
	  double poimin = poi->getMin();
	  double poimax = poi->getMax();
	  //poimin = 0; poimax=10;

	  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
	  calc->SetFixedScan(npoints,poimin,poimax);
 	  calc->SetVerbose(2); 
	  HypoTestInverterResult * r = calc->GetInterval();

	  double upperLimit = r->UpperLimit();

	  std::cout << "The Expected upper limit is: " << r->GetExpectedUpperLimit(0) << std::endl;


/*
//------------ Getting the interval as function of m --------------//
   ifstream in;
   in.open("integral_mass.dat");
   

  vector <double> masses_v;
  vector <double> observed_v;
  vector <double> expected_v;
  vector <double> expected_S1_up_v;
  vector <double> expected_S1_dw_v;
  vector <double> expected_S2_up_v;
  vector <double> expected_S2_dw_v;

  double mass_itr =0.;
  double Nev_exp_th_itr =0.;
  double xsec_modifier = 10.;
  while(mass_itr <1000.){
	in >> mass_itr;
	in >> Nev_exp_th_itr;

 	
	xsec_modifier = Nev_exp_th_itr * 225. * 34.;  //224.6 livedays and 48 kg and 10^-40 cm2 Xsec.

	masses_v.push_back(mass_itr);
	observed_v.push_back( 1e-40 * 1. / xsec_modifier * upperLimit );
	expected_v.push_back( 1e-40  * 1. / xsec_modifier * r->GetExpectedUpperLimit(0) );
	expected_S1_up_v.push_back(1e-40 *1. / xsec_modifier * r->GetExpectedUpperLimit(1));
	expected_S2_up_v.push_back(1e-40 *1. / xsec_modifier * r->GetExpectedUpperLimit(2));
	expected_S2_dw_v.push_back(1e-40 *1. / xsec_modifier * r->GetExpectedUpperLimit(-2));
	expected_S1_dw_v.push_back(1e-40 *1. / xsec_modifier * r->GetExpectedUpperLimit(-1));

        cout << "Mass  " << mass_itr << "  Limit   " <<  1e-40  * 1. / xsec_modifier * r->GetExpectedUpperLimit(0)  << endl;
	
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
double exSigma1_l[n];
double exSigma1_u[n];
double exSigma2_l[n];
double exSigma2_u[n];

for(int k=0; k< n; k++){

	mA[k] = masses_v[k];
	observed[k] = observed_v[k];
	expected[k] = expected_v[k];
	exSigma1_l[k] =expected_v[k] -  expected_S1_dw_v[k] ;
 	exSigma1_u[k] = expected_S1_up_v[k] - expected_v[k];
	exSigma2_l[k] = expected_v[k] - expected_S2_dw_v[k];
	exSigma2_u[k] = expected_S2_up_v[k] - expected_v[k] ;
}

TGraphErrors *obs_limits = new TGraphErrors(n, mA, observed);
TGraphErrors *Exp_limits = new TGraphErrors(n, mA, expected );
TGraphAsymmErrors *Exp_limitsS1 = new TGraphAsymmErrors(n, mA, expected ,xe, xe, exSigma1_l, exSigma1_u );
TGraphAsymmErrors *Exp_limitsS2 = new TGraphAsymmErrors(n, mA, expected ,xe, xe, exSigma2_l, exSigma2_u);

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


//Exp_limitsS2->GetYaxis()->SetTitle("#sigma#timesBR( #phi #rightarrow #tau#tau )  [pb]");
Exp_limitsS2->GetYaxis()->SetTitle("#sigma");

Exp_limitsS2->GetXaxis()->SetTitle("M  [GeV]");


//Exp_limitsS2->GetXaxis()->SetRangeUser(10.,1000.);
//Exp_limitsS2->GetYaxis()->SetRangeUser(1E-38,1E-30);
Exp_limitsS2->GetXaxis()->SetLimits(9.,1000.);
Exp_limitsS2->GetYaxis()->SetRangeUser(1E-38,1E-30);

Exp_limitsS2->Draw("Al3");
Exp_limitsS1->Draw("sameL3");
Exp_limits->Draw("PL");
obs_limits->Draw("PL");


TLegend* lego = new TLegend(0.2,0.9,0.5,0.7);
  lego->SetTextSize(0.033);
  lego->SetFillColor(0);
  lego->SetBorderSize(0);
  lego->AddEntry(obs_limits,"Observed 95\% CLs limit");
  lego->AddEntry(Exp_limits, "Expected 95\% CLs limit");
  lego->AddEntry(Exp_limitsS1,"1 #sigma","f");
  lego->AddEntry(Exp_limitsS2,"2 #sigma","f");
  lego->Draw();


gPad->SetLogy();
gPad->SetLogx();
gPad->RedrawAxis();

myText(0.4,0.86,2,"Test");

*/

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





