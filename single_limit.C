
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TLegend.h"
#include <cstdlib>
#include <cstdio>

/*#include "RooStats/ConfInterval.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/HybridCalculator.h"
#include "RooStats/NumEventsTestStat.h"
#include "RooStats/HypoTestPlot.h"
*/


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
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"
#include "RooStats/ProfileInspector.h"

#include <iostream>     // std::cout, std::ios
#include <fstream>     // std::cout, std::ios
#include <stdio.h>
#include "Math/MinimizerOptions.h"




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
/*
Events in CR data  836
Events in CR Co60  22017
Events in CR AmBe  628
Norm Factor        0.0379707  +-  0.00133794  (rel in % = 3.52362 )
Events in SR Co60  24200
Events in SR AmBe  13083
Events in TOTAL for AmBe  111752
Events in TOTAL for Co60  72157
Events in TOTAL for DATA  2685
TOTAL integrated livetime 225.009 days
*/

for(int step=0; step < 12; step++){

 cout << "STARTING ---------------- " << endl;
  double massPoint =0.;
  double n_bins = 16.;

  double Bkg_norm_factor = 0.0212;//0.0378;
  double err_Bkg_norm_factor = 0.002;
  //double err_Bkg_norm_factor = 0.0001;

  double S_tot 			= 1000000.;//13690;	// Total events in SIGNAL REGIONS for AmBe (--> the Eff_i are the efficiencies relative to SR)
  double B_tot 			= 27030.;//24318;	// Total events in SR for Co60
  double Acceptance_SR_obs	= 0.91;		// Cut Acceptance for Signal
  double err_Acceptance_SR 	= 0.04;
  double N_tot_theory 		= 50.;  	// Total Expected Signal Events for a given mass and Xsec. Set to 10 because gives mu in range [0.,100]

//---- retrieving global variables ----//
  TFile *histos  = TFile::Open("histo_limits_gaudenz_final.root");
  TH1F  *h_ambe_SR;
  if(step ==0 ) {massPoint = 20.;  h_ambe_SR = (TH1F*)histos->Get("hist0");}
  else if(step ==1) {massPoint = 33.;  h_ambe_SR = (TH1F*)histos->Get("hist1");}
  else if(step ==2) {massPoint = 54.; h_ambe_SR = (TH1F*)histos->Get("hist2");}
  else if(step==3)  {massPoint = 90.; h_ambe_SR = (TH1F*)histos->Get("hist3");}
  else if(step ==4) {massPoint = 149.; h_ambe_SR = (TH1F*)histos->Get("hist4");}
  else if(step == 5){ massPoint = 246.; h_ambe_SR = (TH1F*)histos->Get("hist5");}
  else if(step ==6) {massPoint = 406.; h_ambe_SR = (TH1F*)histos->Get("hist6");}
  else if(step ==7) {massPoint = 672.; h_ambe_SR = (TH1F*)histos->Get("hist7");}
  else if(step ==8) {massPoint = 1110.; h_ambe_SR = (TH1F*)histos->Get("hist8");}
  else if(step==9 ) {massPoint = 1835.; h_ambe_SR = (TH1F*)histos->Get("hist9");}
  else if(step ==10){massPoint = 3033.; h_ambe_SR = (TH1F*)histos->Get("hist10");}
  else if(step==11) {massPoint = 5010.; h_ambe_SR = (TH1F*)histos->Get("hist11");}
  else {cout << "not found " << endl; exit(100);}


  TH1F  *h_co60_SR = (TH1F*)histos->Get("Co-hist");
  TH1F  *h_DM_SR = (TH1F*)histos->Get("DM-hist");
  /*TFile *histos  = TFile::Open("h_for_limits.root");
  TH1F  *h_ambe_SR = (TH1F*)histos->Get("ambe_SR");
  TH1F  *h_co60_SR = (TH1F*)histos->Get("co60_SR");
  TH1F  *h_DM_SR = (TH1F*)histos->Get("DM_SR");
*/
//  h_ambe_SR->Rebin(2);
//  h_co60_SR->Rebin(2);
//  h_DM_SR->Rebin(2);

/*h_ambe_SR->Print("all");
h_co60_SR->Print("all");
h_DM_SR->Print("ALL");
h_co60_SR->Scale(Bkg_norm_factor);
h_co60_SR->Print("ALL");
*/

RooWorkspace w("w");
//Poi
RooRealVar mu("mu","POI",1.,0.,5.);
//w.import(mu);

//-- constants
RooRealVar o_Acc("o_Acc","acceptance",0.);
o_Acc.setConstant(true);

RooRealVar M_Acc("M_Acc","acceptance",Acceptance_SR_obs);
M_Acc.setConstant(true);

RooRealVar E_Acc("E_Acc","acceptance",err_Acceptance_SR);
E_Acc.setConstant(true);

RooRealVar E_Acc_T("E_Acc_T","acceptance",1.);
E_Acc_T.setConstant(true);


//w.import(Acc);
//w.import(eff);

RooRealVar n_theory("n_theory","expected events from theory",N_tot_theory);
n_theory.setConstant(true);
//w.import(n_theory);

RooRealVar M_Bkg_norm_factor("M_Bkg_norm_factor","mean value",Bkg_norm_factor);
M_Bkg_norm_factor.setConstant(true);
//w.import(M_Bkg_norm_factor);
RooRealVar E_Bkg_norm_factor("E_Bkg_norm_factor","error",err_Bkg_norm_factor);
E_Bkg_norm_factor.setConstant(true);
//w.import(E_Bkg_norm_factor);


//-- Global Observables  -- dummy values since it is rescaled 
RooRealVar o_bkg_norm("o_bkg_norm","observed value auxiliary measur", 0.);
o_bkg_norm.setConstant(true);
//w.import(o_bkg_norm);
RooRealVar o_err_norm("o_err_norm","observed value auxiliary measur", 1.);
o_err_norm.setConstant(true);


//w.import(o_err_norm);

//-- Nuissance Par.
RooRealVar T_Acc("T_Acc","Signal acceptance T value",0.,-5.,5.);
RooRealVar T_Bkg_norm_factor("T_Bkg_norm_factor","Tvalue for Bkg_norm_factor",0.,-5.,5.);  // Tvalue, centered at 0 between -5sigma and +5 sigma  
//w.import(T_Bkg_norm_factor);

//-- functions for Tvalued nuissance
RooFormulaVar f_Bkg_norm_factor("f_Bkg_norm_factor","T_Bkg_norm_factor*E_Bkg_norm_factor + M_Bkg_norm_factor",RooArgSet(T_Bkg_norm_factor,E_Bkg_norm_factor,M_Bkg_norm_factor));
RooFormulaVar f_Acc("f_Acc","T_Acc*E_Acc + M_Acc",RooArgSet(T_Acc,E_Acc,M_Acc));
//w.import(f_Bkg_norm_factor);
//w.import(s_plus_b);


//w.import(n_obs);


// PDFs
RooGaussian c_bkg_norm("c_bkg_norm","constraint on bkg_norm", o_bkg_norm, T_Bkg_norm_factor, o_err_norm);
RooGaussian c_Acc("c_Acc","constraint on Acc", o_Acc, T_Acc, E_Acc_T);
//w.import(c_bkg_norm);
//w.import(pdf);


RooArgSet observables;
RooArgSet g_obs;
RooArgSet nuissanceP;
RooArgSet product;

  //bin-auxiliary measure variables
  double S_bin = 0.;
  double B_bin =0.;
  double nobs_bin =0.;

for (int bin_itr = 1; bin_itr <= h_ambe_SR->GetNbinsX() ; bin_itr++){


   //retrieve bin info
   S_bin    = h_ambe_SR->GetBinContent(bin_itr);   // Observed events in bin_i for AmBe
   B_bin    = h_co60_SR->GetBinContent(bin_itr);   // Observed events in bin_i for Co60
   nobs_bin = h_DM_SR->GetBinContent(bin_itr); // Observed events in bin_i for DM data

   //cout<<"Efficiency   ---- " <<  S_bin / S_tot  << endl;
   TString bin_name(TString::Itoa(bin_itr, 10));

   //-- stat uncertainty Co
   RooRealVar *M_n_Co = new RooRealVar ("M_n_Co_"+bin_name,"Mean value", B_bin);
   M_n_Co->setConstant(true);
   RooRealVar *E_n_Co = new RooRealVar("E_n_Co_"+bin_name,"Error", sqrt(B_bin));
   E_n_Co->setConstant(true);
   RooRealVar *E_n_Co_T = new RooRealVar("E_n_Co_T_"+bin_name,"Error", 1.);
   E_n_Co_T->setConstant(true);
   RooRealVar *o_n_Co= new RooRealVar("o_n_Co_"+bin_name,"observed value auxiliary measur", 0.);
   o_n_Co->setConstant(true);

   //-- Stat uncertainty AmBe
    

   RooRealVar *M_eff = new RooRealVar("M_eff_"+bin_name,"bin efficiency signal", S_bin / S_tot);
   M_eff->setConstant(true);
   RooRealVar *E_eff = new RooRealVar("E_eff_"+bin_name,"bin Error", sqrt(S_bin) / S_tot);
   E_eff->setConstant(true);
   RooRealVar *E_eff_T = new RooRealVar("E_eff_T_"+bin_name,"bin Error T value", 1.);
   E_eff_T->setConstant(true);
   RooRealVar *o_eff = new RooRealVar("o_eff_"+bin_name,"bin observed", 0.);
   o_eff->setConstant(true);
   

   //-- Observable
   RooRealVar *n_obs = new RooRealVar("n_obs_"+bin_name,"observed events", nobs_bin);
   //RooRealVar *n_obs = new RooRealVar("n_obs_"+bin_name,"observed events", B_bin*Bkg_norm_factor +10./h_ambe_SR->GetNbinsX());

   RooRealVar *T_n_Co = new RooRealVar("T_n_Co_"+bin_name,"T_n_Co_"+bin_name,0.,-5,5.);
   RooRealVar *T_eff = new RooRealVar("T_eff_"+bin_name,"T_eff_"+bin_name,0.,-5,5.);

   RooGaussian *c_eff = new RooGaussian("c_eff_"+bin_name,"constraint eff", *o_eff, *T_eff, *E_eff_T );
   RooGaussian *c_n_Co = new RooGaussian("c_n_Co_"+bin_name,"constraint n_co", *o_n_Co, *T_n_Co, *E_n_Co_T );
   RooFormulaVar *f_n_Co = new RooFormulaVar("f_n_Co_"+bin_name,"Rescaling the T_value", "T_n_Co_"+bin_name+"*E_n_Co_"+bin_name+" + M_n_Co_"+bin_name, RooArgSet(*T_n_Co,*E_n_Co,*M_n_Co));
   RooFormulaVar *f_eff = new RooFormulaVar("f_eff_"+bin_name,"Rescaling the T_value", "T_eff_"+bin_name+"*E_eff_"+bin_name+" + M_eff_"+bin_name, RooArgSet(*T_eff,*E_eff,*M_eff));
   RooFormulaVar *s_plus_b = new RooFormulaVar("s_plus_b_"+bin_name,"n_theory*f_Acc*f_eff_"+bin_name +"*mu + f_n_Co_"+bin_name+"*f_Bkg_norm_factor",RooArgSet(n_theory,*f_eff,f_Acc, mu,*f_n_Co,f_Bkg_norm_factor));

   RooPoisson  *pdf = new RooPoisson("pdf_"+bin_name,"",*n_obs, *s_plus_b);
  
   observables.add(*n_obs);
   g_obs.add(*o_n_Co);
   g_obs.add(*o_eff);
   nuissanceP.add(*T_n_Co);
   nuissanceP.add(*T_eff);

   product.add(*pdf);
   product.add(*c_n_Co);
   product.add(*c_eff);
  

}


histos->Close();
   
   //constraint, stat uncertainty
//   w.factory("Gamma::gamma_"+bin_name+"(N_Co_SR_bin_"+bin_name+",sum::temp_"+bin_name+"(Cobs_"+bin_name+"[1000],1),1,0)") ;
//   w.factory("Poisson:Co_constrain_"+bin_name+"(Cobs_"+bin_name+"[200],N_Co_SR_bin_"+bin_name+")");
   
product.add(c_bkg_norm);
product.add(c_Acc);
RooProdPdf model("model","model",product) ;

//w.import(model);





   g_obs.add(o_bkg_norm);
   g_obs.add(o_Acc);
   nuissanceP.add(T_Bkg_norm_factor);
   nuissanceP.add(T_Acc);

//w.import(observables);
//w.import(g_obs);
//w.import(nuissanceP);

// Building the model
   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(model);
   mc.SetParametersOfInterest(mu);
   //mc.SetPriorPdf(*w.pdf("prior_s"));

// Setting nuissance parameter
   		//RooArgSet nuissanceParameter;  // No Sys
   mc.SetNuisanceParameters(nuissanceP);

//   mc.SetPriorPdf(*w.pdf("prior"));
// need now to set the global observable
  		 //RooArgSet g_obs;  // noSYS
   mc.SetGlobalObservables(g_obs);

   mc.SetObservables(observables);


// this is needed for the hypothesis tests
   mc.SetSnapshot(mu);


// make data set with the number of observed events
RooDataSet data("data","", observables);
data.add(observables);


w.Print();



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
//	  ac.SetQTilde(true);
	  ac.SetPrintLevel(0);  // to suppress print level 

	  

	// create hypotest inverter 
	  // passing the desired calculator 
	  HypoTestInverter *calc = new HypoTestInverter(ac);    // for asymptotic 
	  //HypoTestInverter calc(fc);  // for frequentist


	  calc->SetConfidenceLevel(0.90);
	  calc->UseCLs(true);
	  int npoints = 100;  // number of points to scan
	  //int npoints = 1000;  // number of points to scan
	  // min and max (better to choose smaller intervals)
	  double poimin = poi->getMin();
	  double poimax = poi->getMax();
	  //poimin = 0; poimax=10;

	  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
	  calc->SetFixedScan(npoints,poimin,poimax);
 	  calc->SetVerbose(0); 
	  HypoTestInverterResult * r = calc->GetInterval();

	  double upperLimit = r->UpperLimit();

	  std::cout << "The Expected upper limit is: " << r->GetExpectedUpperLimit(0) << std::endl;
	  std::cout << "The Expected upper +1sigma" << r->GetExpectedUpperLimit(1) << std::endl;
	  std::cout << "The Expected upper -1sigma" << r->GetExpectedUpperLimit(-1) << std::endl;







//------------ Getting the interval as function of m --------------//
   ifstream in;
   in.open("integral_mass_long.dat");
   

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
//  double N_tot_theory = w.var("N_tot_theory")->getValV();

 cout << "STARTING mass loop ---------------- for mass "<< massPoint << endl;
  while(mass_itr < massPoint){  // as soon as you find the Nevent of the mass point you want then stop.

	in >> mass_itr;
	in >> Nev_exp_th_itr;

  }
 cout << "END mass loop ---------------- " << endl;

 in.close();	

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

char a[10] = "";

sprintf(a,"%d",(int)massPoint);
TString name_mass(a);
TFile f("limit_"+name_mass+".root","RECREATE");

obs_limits->Write("obs_limits");
Exp_limits->Write("Exp_limits");
Exp_limitsS1->Write("Exp_limitsS1");
Exp_limitsS2->Write("Exp_limitsS2");

f.Close();
 cout << "END  loop ---------------- " << endl;
}

}
