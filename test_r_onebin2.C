
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TLegend.h"

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

double n_bins = 1.;

double N_bkg_SR = 24318.;
double Bkg_norm_factor = 0.0378;
double err_Bkg_norm_factor = 0.0013;

  double S_tot 			= 13690;	// Total events in SIGNAL REGIONS for AmBe (--> the Eff_i are the efficiencies relative to SR)
  //double S_tot 			= 123552;	// Total events in SIGNAL REGIONS for AmBe (--> the Eff_i are the efficiencies relative to SR)
  double B_tot 			= 24318;	// Total events in SR for Co60
  double Acceptance_SR_obs	= 0.9;		// Cut Acceptance for Signal
  double err_Acceptance_SR 	= 0.05;
  //double err_Acceptance_SR 	= 0.05;
  double N_tot_theory 		= 50.;  	// Total Expected Signal Events for a given mass and Xsec. Set to 10 because gives mu in range [0.,100]



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

RooRealVar eff("eff","bin efficiency signal", 1./n_bins);
eff.setConstant(true);

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

for (int bin_itr = 1; bin_itr <= n_bins ; bin_itr++){


   //retrieve bin info
  /*  Default
   S_bin = 500.;   // Observed events in bin_i for AmBe
   B_bin = 100.;   // Observed events in bin_i for Co60
   nobs_bin = 45.; // Observed events in bin_i for DM data
   */


   TString bin_name(TString::Itoa(bin_itr, 10));

   RooRealVar *M_n_Co = new RooRealVar ("M_n_Co_"+bin_name,"Mean value", N_bkg_SR/n_bins);
   M_n_Co->setConstant(true);
   RooRealVar *E_n_Co = new RooRealVar("E_n_Co_"+bin_name,"Error", sqrt(N_bkg_SR/n_bins));
   E_n_Co->setConstant(true);
   RooRealVar *E_n_Co_T = new RooRealVar("E_n_Co_T","Error", 1.);
   E_n_Co_T->setConstant(true);
   RooRealVar *o_n_Co= new RooRealVar("o_n_Co_"+bin_name,"observed value auxiliary measur", 0.);
   o_n_Co->setConstant(true);


   //-- Observable
   RooRealVar *n_obs = new RooRealVar("n_obs_"+bin_name,"observed events", N_bkg_SR*Bkg_norm_factor/n_bins + 10./n_bins);

   RooRealVar *T_n_Co = new RooRealVar("T_n_Co_"+bin_name,"T_n_Co_"+bin_name,0.,-5,5.);
   RooGaussian *c_n_Co = new RooGaussian("c_n_Co_"+bin_name,"constraint n_co", *o_n_Co, *T_n_Co, *E_n_Co_T );
   RooFormulaVar *f_n_Co = new RooFormulaVar("f_n_Co_"+bin_name,"Rescaling the T_value", "T_n_Co_"+bin_name+"*E_n_Co_"+bin_name+" + M_n_Co_"+bin_name, RooArgSet(*T_n_Co,*E_n_Co,*M_n_Co));
   RooFormulaVar *s_plus_b = new RooFormulaVar("s_plus_b_"+bin_name,"n_theory*eff*f_Acc*mu + f_n_Co_"+bin_name+"*f_Bkg_norm_factor",RooArgSet(n_theory,eff,f_Acc, mu,*f_n_Co,f_Bkg_norm_factor));

   RooPoisson  *pdf = new RooPoisson("pdf_"+bin_name,"",*n_obs, *s_plus_b);
  
   observables.add(*n_obs);
   g_obs.add(*o_n_Co);
   nuissanceP.add(*T_n_Co);

   product.add(*pdf);
   product.add(*c_n_Co);
  

}


   
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


  //////////////////////////////////////////////
  // now use the profile inspector
  ProfileInspector p;
  TList* list = p.GetListOfProfilePlots(data,&mc);
  
  // now make plots
  TCanvas* c1 = new TCanvas("c1","ProfileInspectorDemo",800,600);
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

w.Print();

/*
RooPlot *framei2 = mu.frame();
model.plotOn(framei2);
framei2->Draw();

new TCanvas();
RooPlot *frame = T_Bkg_norm_factor.frame();
model.plotOn(frame);
frame->Draw();
*/


/*
w.factory("Acc[0.9]");    
w.factory("eff[1.]");
w.factory("mu[1.,-1.,50.]");
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
w.var("err_Norm_factor")->setVal(0.0013); 		/// Set error on norm factor   ////SEMBRA CHE IL PROBLEMA SIA PROPRIO NELL'ERRORE DI QUESTO

w.var("N_tot_theory")->setConstant(true); 
w.var("Bkg_norm_factor_obs")->setConstant(true); 
w.var("err_Norm_factor")->setConstant(true); 

w.var("Bkg_norm_factor")->setMax(w.var("Bkg_norm_factor_obs")->getValV() + w.var("err_Norm_factor")->getValV()*5.); //set to 10 sigma.
w.var("Bkg_norm_factor")->setMin(w.var("Bkg_norm_factor_obs")->getValV() - w.var("err_Norm_factor")->getValV()*5.); //set to 10 sigma.
w.var("Bkg_norm_factor")->setVal(w.var("Bkg_norm_factor_obs")->getValV());// Set to observed value!!!

for(int bin_itr = 1; bin_itr <=n_bins; bin_itr++){

   TString bin_name(TString::Itoa(bin_itr, 10));
   w.factory("prod:N_s_"+bin_name+"(N_tot_theory, Acc, eff, mu)");
   w.factory("prod:N_b_"+bin_name+"(N_Co_SR_bin_"+bin_name+"[240,120,360], Bkg_norm_factor)");
   w.factory("sum:nexp_"+bin_name+"(N_s_"+bin_name+", N_b_"+bin_name+")");

   w.factory("Poisson:pdf_"+bin_name+"(nobs_"+bin_name+"[10],nexp_"+bin_name+")");

   //constraint, stat uncertainty
//   w.factory("Gamma::gamma_"+bin_name+"(N_Co_SR_bin_"+bin_name+",sum::temp_"+bin_name+"(Cobs_"+bin_name+"[1000],1),1,0)") ;
//   w.factory("Poisson:Co_constrain_"+bin_name+"(Cobs_"+bin_name+"[200],N_Co_SR_bin_"+bin_name+")");
   
   w.factory("err_"+bin_name+"[500]");
   w.var("err_"+bin_name)->setVal(sqrt(N_bkg_SR/n_bins));
   w.var("err_"+bin_name)->setConstant(true);
	
   w.factory("Gaussian:Co_constrain_"+bin_name+"(Cobs_"+bin_name+"[200],N_Co_SR_bin_"+bin_name+",err_"+bin_name+")");


   w.var("nobs_"+bin_name)->setVal(N_bkg_SR*0.0378/n_bins + 1.);
   w.var("N_Co_SR_bin_"+bin_name)->setVal(N_bkg_SR/n_bins);
   w.var("N_Co_SR_bin_"+bin_name)->setMax(N_bkg_SR/n_bins + sqrt(N_bkg_SR/n_bins)*5.);
   w.var("N_Co_SR_bin_"+bin_name)->setMin(N_bkg_SR/n_bins - sqrt(N_bkg_SR/n_bins)*5.);
   w.var("Cobs_"+bin_name)->setVal(N_bkg_SR/n_bins);
   w.var("Cobs_"+bin_name)->setConstant(true); 

  //w.var("N_Co_SR_bin_"+bin_name)->setConstant(true);  

  // total_model.Append("pdf_"+bin_name+",gamma_"+bin_name+",");
   total_model.Append("pdf_"+bin_name+",Co_constrain_"+bin_name+",");

   observables.add(*w.var("nobs_"+bin_name));

   g_obs.add(*w.var("Cobs_"+bin_name));
   nuissanceP.add(*w.var("N_Co_SR_bin_"+bin_name));

} 

//------------priors   
//  w.factory("Uniform::prior_nuis({Bkg_norm_factor,N_Co_SR_bin_1})");
//  w.factory("Uniform::prior_nuis({Bkg_norm_factor,N_Co_SR_bin_1,N_Co_SR_bin_2})");
//  w.factory("Uniform::prior_poi({mu})");
//  w.factory("PROD::prior(prior_poi,prior_nuis)"); 

total_model.Append("costraint_bkg_norm)");
w.factory(total_model);

w.Print();
*/


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

/*
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

   //w.writeToFile("test.root", true);




w.Print();*/

/*
// example use of Feldman-Cousins
  FeldmanCousins fc(data, mc); 
  fc.SetConfidenceLevel( 0.9);
  fc.SetNBins(100); // number of points to test per parameter
  fc.UseAdaptiveSampling(true); // make it go faster

  // Here, we consider only ensembles with 100 events
  // The PDF could be extended and this could be removed
  fc.FluctuateNumDataEntries(false); 

  // Proof
  //  ProofConfig pc(*wspace, 4, "workers=4", kFALSE);    // proof-lite
  //ProofConfig pc(w, 8, "localhost");    // proof cluster at "localhost"
  //  ToyMCSampler* toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler();
  //  toymcsampler->SetProofConfig(&pc);     // enable proof

  PointSetInterval* interval = (PointSetInterval*) fc.GetInterval();

 // example use profile likelihood calculator
  ProfileLikelihoodCalculator plc(data, mc);
  plc.SetConfidenceLevel( 0.9);
  LikelihoodInterval* plInt = plc.GetInterval();


  RooRealVar* mu = w.var("mu");

 cout << "plc interval is [" << 
    plInt->LowerLimit(*mu) << ", " << 
    plInt->UpperLimit(*mu) << "]" << endl;

std::cout << "fc interval is ["<<  
    interval->LowerLimit(*mu) << " , "  <<
    interval->UpperLimit(*mu) << "]" << endl;
*/
//data.Print();
/*
  //////////////////////////////////////////////
  // now use the profile inspector
  ProfileInspector p;
  TList* list = p.GetListOfProfilePlots(data,&mc);
  
  // now make plots
  TCanvas* c1 = new TCanvas("c1","ProfileInspectorDemo",800,600);
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



*/

//////////////// hybrid calculator ----//
/*
 HybridCalculator hc1(data, *sbModel, *bModel);
  ToyMCSampler *toymcs1 = (ToyMCSampler*)hc1.GetTestStatSampler();
  //  toymcs1->SetNEventsPerToy(1); // because the model is in number counting form
  toymcs1->SetTestStatistic(&binCount); // set the test statistic
  //  toymcs1->SetGenerateBinned();
  hc1.SetToys(30000,1000); 
//  hc1.ForcePriorNuisanceAlt(*w->pdf("py"));
//  hc1.ForcePriorNuisanceNull(*w->pdf("py"));
  // if you wanted to use the ad hoc Gaussian prior instead
  //  hc1.ForcePriorNuisanceAlt(*w->pdf("gauss_prior"));
  //  hc1.ForcePriorNuisanceNull(*w->pdf("gauss_prior"));
  // if you wanted to use the ad hoc log-normal prior instead
  //  hc1.ForcePriorNuisanceAlt(*w->pdf("lognorm_prior"));
  //  hc1.ForcePriorNuisanceNull(*w->pdf("lognorm_prior"));

  // enable proof
  // proof not enabled for this test statistic
  //  if(pc) toymcs1->SetProofConfig(pc);     

  // these lines save current msg level and then kill any messages below ERROR
//  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  // Get the result
  HypoTestResult *r1 = hc1.GetHypoTest();
  //RooMsgService::instance().setGlobalKillBelow(msglevel); // set it back
  cout << "-----------------------------------------"<<endl;
  cout << "Part 4" << endl;
  r1->Print();
 // t.Stop();  t.Print(); t.Reset(); t.Start();

 // c->cd(2);
  HypoTestPlot *p1 = new HypoTestPlot(*r1,30); // 30 bins, TS is discrete
t  p1->Draw();
*/

/*
//----------------   Limit testing   -----------------///
      TCanvas* c2 = new TCanvas("c2","NuissanceTest",800,600);
      bool IsConditionnal = false;
	
      RooAbsPdf *pdftmp = mc.GetPdf();
      ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
      RooFitResult *fitres = pdftmp->fitTo( data , Save() , Minos(false) );
      cout << endl;
      cout << endl;
      if (IsConditionnal) cout << "Conditionnal fit : mu is fixed at " << endl; //<< mu << endl;
      else                cout << "Unconditionnal fit : mu is fitted" << endl;
   //   double muhat = firstPOI->getVal();
   //   firstPOI->setConstant(kFALSE);
  
    //--- correlation matrix 
    TH2D *h2Dcorrelation = (TH2D*) fitres->correlationHist();
    h2Dcorrelation->Draw("colz");
  

   // Plotting the likelihood projection in each NP direction
    RooRealVar* var = NULL;
    RooRealVar * poi_test =  (RooRealVar*)(mc.GetParametersOfInterest()->first());
    RooAbsReal* nll = pdftmp->createNLL(data);
    TIterator* it3 = mc.GetNuisanceParameters()->createIterator();
    while( (var = (RooRealVar*) it3->Next()) ){
      TString vname=var->GetName();
      RooPlot* frame2 = var->frame(Title("-log(L) vs "+vname)) ;
      nll->plotOn(frame2,LineColor(kRed),ShiftToZero()) ;
      frame2->GetYaxis()->SetRangeUser(0.0,5.0);
      frame2->GetYaxis()->SetTitle("#Delta [-2Log(L)]");
      new TCanvas( "NLLscan_"+vname );
      frame2->Draw();
   }
     RooPlot* frame2 = poi_test->frame(Title("-log(L) vs POI")) ;
     nll->plotOn(frame2,LineColor(kRed),ShiftToZero()) ;frame2->GetYaxis()->SetRangeUser(0.0,5.0);
      frame2->GetYaxis()->SetTitle("#Delta [-2Log(L)]");
      new TCanvas( "NLLscan_mu" );
      frame2->Draw();

     
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


//------------------Limit calculation for N_th event expected = 1


	  AsymptoticCalculator  ac(data, *bModel, *sbModel);
	  //ac.SetOneSidedDiscovery(true);  // for one-side discovery test
	  ac.SetOneSided(true);  // for one-side tests (limits)
//	  ac.SetQTilde(true);
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
	  std::cout << "The Expected upper +1sigma" << r->GetExpectedUpperLimit(1) << std::endl;
	  std::cout << "The Expected upper -1sigma" << r->GetExpectedUpperLimit(-1) << std::endl;









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





