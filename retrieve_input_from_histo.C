
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


void retrieve_input_from_histo(RooWorkspace &w){



// Variable definition
   /*
  */

/*********  GLOBAL VARIABLES   ********************/

/*   Default ---
  double Bkg_norm_factor_obs  	= 0.5 ;		// Normalizzation factor for background (N_Data_CR / N_Co_CR)
  double err_Norm_factor_obs  	= 0.05;
  double S_tot 			= 1000.;	// Total events in SR for AmBe
  double B_tot 			= 10000.;	// Total events in SR for Co60
  double Acceptance_SR_obs	= 0.9;		// Cut Acceptance for Signal
  double err_Acceptance_SR 	= 0.01;
  double N_tot_theory 		= 5.;  	   	// Total Expected Signal Events for a given mass and Xsec.
   */

  double Bkg_norm_factor_obs  	= 0.0378 ;	// Normalizzation factor for background (N_Data_CR / N_Co_CR)
  //double err_Norm_factor_obs  	= 0.00013;
  double err_Norm_factor_obs  	= 0.0013;
  double S_tot 			= 13690;	// Total events in SIGNAL REGIONS for AmBe (--> the Eff_i are the efficiencies relative to SR)
  //double S_tot 			= 123552;	// Total events in SIGNAL REGIONS for AmBe (--> the Eff_i are the efficiencies relative to SR)
  double B_tot 			= 24318;	// Total events in SR for Co60
  double Acceptance_SR_obs	= 0.9;		// Cut Acceptance for Signal
  double err_Acceptance_SR 	= 0.005;
  //double err_Acceptance_SR 	= 0.05;
  double N_tot_theory 		= 1.;  	// Total Expected Signal Events for a given mass and Xsec. Set to 10 because gives mu in range [0.,100]
/**************************************************/  



//---- retrieving global variables ----//
  TFile *histos  = TFile::Open("h_for_limits.root");
  TH1D  *h_ambe_SR = (TH1D*)histos->Get("ambe_SR_flat");
  TH1D  *h_co60_SR = (TH1D*)histos->Get("co60_SR_flat");
  TH1D  *h_DM_SR = (TH1D*)histos->Get("DM_SR_flat");

  h_ambe_SR->Rebin(10);
  h_co60_SR->Rebin(10);
  h_DM_SR->Rebin(10);



//----- Model for Global Variables  -----//
   w.factory("N_tot_theory[10]");
   w.factory("Acceptance_SR[0.9,0.0,1]");
   w.factory("mu[1.,0.,100]");
   w.factory("Bkg_norm_factor[0.7, 0.0,10.0]");
   w.factory("S_tot[1000.0]");


   w.factory("Gaussian:costraint_bkg_norm(Bkg_norm_factor_obs[0,10], Bkg_norm_factor, err_Norm_factor[1])"); // approximation (should be kind of complicated cauchy....)
   w.factory("Gaussian:costarint_sig_acc(Acceptance_SR_obs[0,1], Acceptance_SR, err_Acceptance_SR[0.01])");// this has the problem of boundaries acceptance cannot be larger than 1.  --->Lognormal???

   RooArgSet nuissanceParameter (*w.var("Bkg_norm_factor"),*w.var("Acceptance_SR"));
   RooArgSet g_obs(*w.var("Bkg_norm_factor_obs"), *w.var("Acceptance_SR_obs"));
   RooArgSet obs;

//----- Setting Value to Global Variables ------//
   w.var("N_tot_theory")->setVal(N_tot_theory); 		// Total expected Signal events
   w.var("Bkg_norm_factor_obs")->setVal(Bkg_norm_factor_obs); 	// Ratio between data in CR and Co in CR
   w.var("err_Norm_factor")->setVal(err_Norm_factor_obs); 	// Set error on norm factor
   w.var("Acceptance_SR_obs")->setVal(Acceptance_SR_obs); 
   w.var("err_Acceptance_SR")->setVal(err_Acceptance_SR);
   w.var("S_tot")->setVal(S_tot); 				// Total events of AmBe in SR 

   w.var("N_tot_theory")->setConstant(true); 
   w.var("S_tot")->setConstant(true); 
   w.var("Bkg_norm_factor_obs")->setConstant(true); 
   w.var("err_Norm_factor")->setConstant(true); 
   w.var("Acceptance_SR_obs")->setConstant(true); 
   w.var("err_Acceptance_SR")->setConstant(true); 

//-------- Setting Value of NP to Observed!!   ---------------//
   w.var("Acceptance_SR")->setVal(w.var("Acceptance_SR_obs")->getValV());
   w.var("Acceptance_SR")->setMax(w.var("Acceptance_SR_obs")->getValV() + w.var("err_Acceptance_SR")->getValV() * 5.);// Set to 10 sigma.
   w.var("Bkg_norm_factor")->setMax(w.var("Bkg_norm_factor_obs")->getValV() + w.var("err_Norm_factor")->getValV()*5.); //set to 10 sigma.
   w.var("Bkg_norm_factor")->setVal(w.var("Bkg_norm_factor_obs")->getValV());// Set to observed value!!!
  
  
  TString total_model = "PROD:model(";

  //bin-auxiliary measure variables
  double S_bin = 0.;
  double B_bin =0.;
  double nobs_bin =0.;
  
 /*****************  LOOP OVER  BINS   *******************/
  for (int bin_itr = 1; bin_itr <= h_ambe_SR->GetNbinsX() ; bin_itr++){


   //retrieve bin info
  /*  Default
   S_bin = 500.;   // Observed events in bin_i for AmBe
   B_bin = 100.;   // Observed events in bin_i for Co60
   nobs_bin = 45.; // Observed events in bin_i for DM data
   */
   S_bin    = h_ambe_SR->GetBinContent(bin_itr);   // Observed events in bin_i for AmBe
   B_bin    = h_co60_SR->GetBinContent(bin_itr);   // Observed events in bin_i for Co60
   nobs_bin = h_DM_SR->GetBinContent(bin_itr); // Observed events in bin_i for DM data

   TString bin_name(TString::Itoa(bin_itr, 10));

   w.factory("prod:N_s_"+bin_name+"(N_tot_theory, Acceptance_SR, Eff_AmBe_bin_"+bin_name+"[0.01,0.0,1], mu)");
   w.factory("prod:N_b_"+bin_name+"(N_Co_SR_bin_"+bin_name+"[100.0,0.0,50000.0], Bkg_norm_factor)"); 
   w.factory("sum:nexp_"+bin_name+"(N_s_"+bin_name+", N_b_"+bin_name+")");
// Poisson of (n | s+b)
   w.factory("Poisson:pdf_"+bin_name+"(nobs_"+bin_name+"[500,0.0,5000.0],nexp_"+bin_name+")");


// Poisson constraint --- Stat uncertainty on the bin content
   w.factory("prod:S_"+bin_name+"_exp(S_tot, Eff_AmBe_bin_"+bin_name+")");   // put 1000. to right AmBe amount
   w.factory("Poisson:AmBe_bin_"+bin_name+"(S_"+bin_name+"[500,0.0,50000.0],S_"+bin_name+"_exp)");  // change 50 for multiple bin, put 1000. to right AmBe amount
   w.factory("Poisson:Co_SR_bin_"+bin_name+"(B_"+bin_name+"[500,0.0,50000.0], N_Co_SR_bin_"+bin_name+")");  //change 100 for multiple bin, put 1000. to right Co tot amount in SR

   w.factory("PROD:model_"+bin_name+"(pdf_"+bin_name+",AmBe_bin_"+bin_name+",Co_SR_bin_"+bin_name+")");
// Set input variables -- Global OBSERVABLESerr_Norm_factor
   w.var("S_"+bin_name)->setMax(w.var("S_tot")->getValV());  // the range of the poisson needs to be defined
   w.var("S_"+bin_name)->setVal(S_bin); 	  		//  Events in bin i for AmBe

   w.var("B_"+bin_name)->setMax(B_tot);  // set to TOT Co in SR
   w.var("B_"+bin_name)->setVal(B_bin);    // Set observed event in SR Co
 
   w.var("S_"+bin_name)->setConstant(true); 
   w.var("B_"+bin_name)->setConstant(true); 


// Set range and value for nuissance parameters
   w.var("N_Co_SR_bin_"+bin_name)->setMax(B_bin + sqrt(B_bin)*5.); //set to 10 sigma.
   w.var("N_Co_SR_bin_"+bin_name)->setVal(B_bin);  		// set to observed value!!
   w.var("Eff_AmBe_bin_"+bin_name)->setVal(S_bin/ S_tot ); // Set to observed value!!!
   w.var("Eff_AmBe_bin_"+bin_name)->setMax( (S_bin + sqrt(S_bin)*5.) / S_tot ); // Set to observed value!!!
   if(w.var("Eff_AmBe_bin_"+bin_name)->getMax() > 1.) cout << "WARNING!!! Eff_AmBe_bin_"+bin_name+" efficiency > 1. " << endl;

//Set Observed events!
   w.var("nobs_"+bin_name)->setVal(nobs_bin);

   nuissanceParameter.add(*w.var("N_Co_SR_bin_"+bin_name));
   nuissanceParameter.add(*w.var("Eff_AmBe_bin_"+bin_name));
   g_obs.add(*w.var("S_"+bin_name));
   g_obs.add(*w.var("B_"+bin_name));

   obs.add(*w.var("nobs_"+bin_name));

   total_model.Append("model_"+bin_name+", ");

} 
/************************************************************************/

total_model.Append("costraint_bkg_norm, costarint_sig_acc)");

w.factory(total_model);

w.defineSet("nuissance_parameter",nuissanceParameter);
w.defineSet("observables",obs);
w.defineSet("g_observables",g_obs);


histos->Close();

  // RooArgSet set = w.allPdfs();
//   set.Print();
// Set input Observables
//  w.var("nobs_"+bin_name+"")->setVal(40);

/*   // various printing
   w.pdf("AmBe_bin_i")->Print("all");

// Building the model
   ModelConfig mc("ModelConfig",&w);
   mc.SetPdf(*w.pdf("model"));
   mc.SetParametersOfInterest(*w.var("mu"));
   //mc.SetPriorPdf(*w.pdf("prior_s"));

// Setting nuissance parameter
   RooArgSet nuissanceParameter (*w.var("Bkg_norm_factor"),*w.var("Acceptance_SR"));
   nuissanceParameter.add(*w.var("N_Co_SR_bin_i"));
   nuissanceParameter.add(*w.var("Eff_AmBe_bin_i"));
   mc.SetNuisanceParameters(nuissanceParameter);

// need now to set the global observable
   RooArgSet g_obs(*w.var("S_i"), *w.var("B_i"));
//   g_obs.add(*w.var("N_tot_theory"));
//   g_obs.add(*w.var("S_tot"));
   g_obs.add(*w.var("Bkg_norm_factor_obs"));
//   g_obs.add(*w.var("err_Norm_factor"));
   g_obs.add(*w.var("Acceptance_SR_obs"));
//   g_obs.add(*w.var("err_Acceptance_SR"));
   mc.SetGlobalObservables(g_obs);

   RooArgSet observables(*w.var("nobs_i"));
   //observables.add(variousbins...);
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

   w.writeToFile("CountingModel.root", true);
*/


}
