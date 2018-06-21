
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "TLegend.h"
#include <cstdlib>
#include <cstdio>


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




void limit_inelastic_stable(TString massPoint) {

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
// CR double sided contains:
 CO60 low 14495  DM low  526  tau = 0.0362884
 CO60 high 10958  DM high  391  tau = 0.0356817
 TAU  = 0.035985  +-  4.43861%  OR 0.00159724

 DM 		in SR 764
 Bkg Expected   in SR 800
 BKG total      in SR 22225

 TOTAL integrated livetime 225.009 days and 34Kg
*/


 cout << "----------- STARTING Inelastic Limits  on mass point:    " << massPoint + " GeV" <<  endl;

  double Bkg_norm_factor        = 0.034;
  double err_Bkg_norm_factor    = 0.002;

  double err_Acceptance_SR 	= 0.07;         // 7% acc uncertainty  -- including SYS variatons
  double N_tot_theory 		= 50.;  	// Total Expected Signal Events, usefull to have always the histo scale to this number. 

//---- retrieving global variables ----//
  TFile *histos  = TFile::Open("DATA/merged_signal_bkg.root");

  TH1D  *h_signal_SR = (TH1D*)histos->Get("PL_SR_s1s2_mass_"+massPoint+"_sys");
  TH1F  *h_co60_SR   = (TH1F*)histos->Get("co60_PL_SR_reb");
  TH1F  *h_DM_SR     = (TH1F*)histos->Get("dm_PL_SR_reb");
/*
  h_signal_SR->Rebin();
  h_co60_SR->Rebin();
  h_DM_SR->Rebin();
*/

  //Scale the signal histo.
  double Signal_yeld_10m40_cm2 = h_signal_SR->Integral();  // N signal events per Day/Kg
  h_signal_SR->Scale(N_tot_theory / Signal_yeld_10m40_cm2);
  

/*h_signal_SR->Print("all");
h_co60_SR->Print("all");
h_DM_SR->Print("ALL");
h_co60_SR->Scale(Bkg_norm_factor);
h_co60_SR->Print("ALL");
*/

RooWorkspace w("w");
//Poi
RooRealVar mu("mu","POI",1.,0.,5.);

//-- constants
/* //not in use anymore
RooRealVar o_Acc("o_Acc","acceptance",0.);                              
o_Acc.setConstant(true);

RooRealVar M_Acc("M_Acc","acceptance", 1.);            
M_Acc.setConstant(true);

//RooRealVar E_Acc("E_Acc","acceptance", err_Acceptance_SR);
E_Acc.setConstant(true);                  //not in use anymore

//RooRealVar E_Acc_T("E_Acc_T","acceptance",1.);
E_Acc_T.setConstant(true);    not in use anymore
*/


RooRealVar M_Bkg_norm_factor("M_Bkg_norm_factor","mean value",Bkg_norm_factor);
M_Bkg_norm_factor.setConstant(true);

RooRealVar E_Bkg_norm_factor("E_Bkg_norm_factor","error",err_Bkg_norm_factor);
E_Bkg_norm_factor.setConstant(true);


//-- Global Observables  -- dummy values since it is rescaled 
RooRealVar o_bkg_norm("o_bkg_norm","observed value auxiliary measur", 0.);
o_bkg_norm.setConstant(true);

RooRealVar o_err_norm("o_err_norm","observed value auxiliary measur", 1.);
o_err_norm.setConstant(true);



//-- Nuissance Par.
//RooRealVar T_Acc("T_Acc","Signal acceptance T value",0.,-5.,5.); //not in use anymore
RooRealVar T_Bkg_norm_factor("T_Bkg_norm_factor","Tvalue for Bkg_norm_factor",0.,-5.,5.);  // Tvalue, centered at 0 between -5sigma and +5 sigma  

//-- functions for Tvalued nuissance
RooFormulaVar f_Bkg_norm_factor("f_Bkg_norm_factor","Bkg norm Factor", "T_Bkg_norm_factor*E_Bkg_norm_factor + M_Bkg_norm_factor",RooArgSet(T_Bkg_norm_factor,E_Bkg_norm_factor,M_Bkg_norm_factor));
//RooFormulaVar f_Acc("f_Acc","Acceptance ","T_Acc*E_Acc + M_Acc",RooArgSet(T_Acc,E_Acc,M_Acc)); //not in use anymore


// PDFs
RooGaussian c_bkg_norm("c_bkg_norm","constraint on bkg_norm", o_bkg_norm, T_Bkg_norm_factor, o_err_norm);
//RooGaussian c_Acc("c_Acc","constraint on Acc", o_Acc, T_Acc, E_Acc_T);   // not in use anymore


RooArgSet observables;
RooArgSet g_obs;
RooArgSet nuissanceP;
RooArgSet product;
//RooArgSet n_S_list;

double    n_S_TOT = 0.;
TString   f_n_S_TOT_str = "";

for (int bin_itr = 1; bin_itr <= h_signal_SR->GetNbinsX() ; bin_itr++){


  //bin-auxiliary measure variables
   double S_bin       = h_signal_SR->GetBinContent(bin_itr);   // Observed events in bin_i for signal model scaled to Ntheory events
   double S_bin_err   = h_signal_SR->GetBinError(bin_itr);     // Error on bin_i correctly scaled
   double B_bin       = h_co60_SR->GetBinContent(bin_itr);     // Observed events in bin_i for Co60
   double B_bin_err   = sqrt( pow(h_co60_SR->GetBinError(bin_itr),2.) + pow(B_bin * 0.04, 2.) ); // stat unc. summed in quadrature with 4% sys uncertainty per bin (coming from Co60 Th232 discrepancy)
   double nobs_bin    = h_DM_SR->GetBinContent(bin_itr);       // Observed events in bin_i for DM data


   n_S_TOT +=  S_bin;   // integral of signal for sys shape constraint

   //cout<<"Efficiency   ---- " <<  S_bin / S_tot  << endl;
   TString bin_name(TString::Itoa(bin_itr, 10));

   //-- stat uncertainty Co
   RooRealVar *M_n_Co = new RooRealVar ("M_n_Co_"+bin_name,"Mean value", B_bin);
   M_n_Co->setConstant(true);

  //if(bin_itr == 9 || bin_itr == 8) B_bin_err = B_bin_err * 2.;

   RooRealVar *E_n_Co = new RooRealVar("E_n_Co_"+bin_name,"Error", B_bin_err );
   E_n_Co->setConstant(true);
   RooRealVar *E_n_Co_T = new RooRealVar("E_n_Co_T_"+bin_name,"Error", 1.);
   E_n_Co_T->setConstant(true);
   RooRealVar *o_n_Co= new RooRealVar("o_n_Co_"+bin_name,"observed value auxiliary measur", 0.);
   o_n_Co->setConstant(true);
   RooRealVar *T_n_Co = new RooRealVar("T_n_Co_"+bin_name,"T_n_Co_"+bin_name,0.,-5,5.);

   //-- Stat uncertainty signal
   RooRealVar *M_n_S = new RooRealVar("M_n_S_"+bin_name,"bin mean val signal", S_bin );  
   M_n_S->setConstant(true);
   RooRealVar *E_n_S = new RooRealVar("E_n_S_"+bin_name,"bin Error signal", S_bin_err );
   E_n_S->setConstant(true);
   RooRealVar *E_n_S_T = new RooRealVar("E_n_S_T_"+bin_name,"bin Error T value", 1.);
   E_n_S_T->setConstant(true);
   RooRealVar *o_n_S = new RooRealVar("o_n_S_"+bin_name,"bin observed", 0.);
   o_n_S->setConstant(true);
   RooRealVar *T_n_S = new RooRealVar("T_n_S_"+bin_name,"T_n_S_"+bin_name,0.,-5,5.);
   

   //-- Observable
   RooRealVar *n_obs = new RooRealVar("n_obs_"+bin_name,"observed events", nobs_bin);
   //RooRealVar *n_obs = new RooRealVar("n_obs_"+bin_name,"observed events", B_bin*Bkg_norm_factor +10./h_signal_SR->GetNbinsX());

   //Stat constraint terms
   RooGaussian *c_n_S = new RooGaussian("c_n_S_"+bin_name,"constraint stat signal", *o_n_S, *T_n_S, *E_n_S_T );
   RooGaussian *c_n_Co = new RooGaussian("c_n_Co_"+bin_name,"constraint n_co", *o_n_Co, *T_n_Co, *E_n_Co_T );

   //Building the N_s + N_b
   RooFormulaVar *f_n_Co = new RooFormulaVar("f_n_Co_"+bin_name,"Rescaling the T_value", "T_n_Co_"+bin_name+"*E_n_Co_"+bin_name+" + M_n_Co_"+bin_name, RooArgSet(*T_n_Co,*E_n_Co,*M_n_Co));
   RooFormulaVar *f_n_S = new RooFormulaVar("f_n_S_"+bin_name,"Rescaling the T_value", "T_n_S_"+bin_name+"*E_n_S_"+bin_name+" + M_n_S_"+bin_name, RooArgSet(*T_n_S,*E_n_S,*M_n_S));
   RooFormulaVar *s_plus_b = new RooFormulaVar("s_plus_b_"+bin_name,"S plus B", "mu*f_n_S_"+bin_name +" + f_n_Co_"+bin_name+"*f_Bkg_norm_factor",RooArgSet(mu,*f_n_S,*f_n_Co,f_Bkg_norm_factor));

   RooPoisson  *pdf = new RooPoisson("pdf_"+bin_name,"",*n_obs, *s_plus_b);
  
   observables.add(*n_obs);
   g_obs.add(*o_n_Co);
   g_obs.add(*o_n_S);
   nuissanceP.add(*T_n_Co);
   nuissanceP.add(*T_n_S);

   product.add(*pdf);
   product.add(*c_n_Co);
   product.add(*c_n_S);
  
/*   if(bin_itr < h_signal_SR->GetNbinsX()) f_n_S_TOT_str.Append("f_n_S_"+bin_name+"+");
   else f_n_S_TOT_str.Append("f_n_S_"+bin_name);
 
   n_S_list.add(*f_n_S);
  */ 	
}


/*RooFormulaVar f_n_S_TOT("f_n_S_TOT","f_n_S_TOT","("+f_n_S_TOT_str+"-50.)/50.",n_S_list); 

RooRealVar M_n_S_TOT("M_n_S_TOT","mean value", 0.);
//RooRealVar M_n_S_TOT("M_n_S_TOT","mean value",n_S_TOT);
M_n_S_TOT.setConstant(true);

RooRealVar E_n_S_TOT("E_n_S_TOT","error", err_Acceptance_SR );
//RooRealVar E_n_S_TOT("E_n_S_TOT","error", n_S_TOT * err_Acceptance_SR );
E_n_S_TOT.setConstant(true);

RooGaussian   c_n_S_TOT("c_n_S_TOT","constraint on N signal",M_n_S_TOT, f_n_S_TOT, E_n_S_TOT); 

g_obs.add(M_n_S_TOT);

product.add(c_n_S_TOT);
*/

histos->Close();
   
   
product.add(c_bkg_norm);
//product.add(c_Acc);        // not in use anymore
RooProdPdf model("model","model",product) ;






   g_obs.add(o_bkg_norm);
   nuissanceP.add(T_Bkg_norm_factor);
//   g_obs.add(o_Acc);         
//   nuissanceP.add(T_Acc);   // not in use anymore



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


//------------------Limit calculation for N_th event expected = 50


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





 // limit writing

  TGraphErrors *obs_limits = new TGraphErrors(1);
  TGraphErrors *Exp_limits = new TGraphErrors(1);
  TGraphAsymmErrors *Exp_limitsS1 = new TGraphAsymmErrors(1);
  TGraphAsymmErrors *Exp_limitsS2 = new TGraphAsymmErrors(1);


  double xsec_modifier = Signal_yeld_10m40_cm2 * 225.009 * 34.;  //225.009 livedays and 34 kg and 10^-40 cm2 Xsec.

  double obs_limit_val = 1.e-40  * N_tot_theory / xsec_modifier * upperLimit;
  double exp_limit_val = 1.e-40  * N_tot_theory / xsec_modifier * r->GetExpectedUpperLimit(0) ; 
  double exp_limit_low_1s = 1.e-40  * N_tot_theory / xsec_modifier * ( r->GetExpectedUpperLimit(0) - r->GetExpectedUpperLimit(-1)); 
  double exp_limit_low_2s = 1.e-40  * N_tot_theory / xsec_modifier * ( r->GetExpectedUpperLimit(0) - r->GetExpectedUpperLimit(-2)); 
  double exp_limit_high_1s = 1.e-40  * N_tot_theory / xsec_modifier * ( r->GetExpectedUpperLimit(1) - r->GetExpectedUpperLimit(0) ); 
  double exp_limit_high_2s = 1.e-40  * N_tot_theory / xsec_modifier * ( r->GetExpectedUpperLimit(2) - r->GetExpectedUpperLimit(0) ); 

  double mass = massPoint.Atof() ;
  obs_limits->SetPoint(0, mass, obs_limit_val );
  Exp_limits->SetPoint(0, mass, exp_limit_val );
  Exp_limitsS1->SetPoint(0, mass, exp_limit_val);
  Exp_limitsS1->SetPointError(0, 0.,0., exp_limit_low_1s, exp_limit_high_1s);
  Exp_limitsS2->SetPoint(0, mass, exp_limit_val);
  Exp_limitsS2->SetPointError(0, 0.,0., exp_limit_low_2s, exp_limit_high_2s);
 
  cout << "Expected median limit for mass " << massPoint << "GeV  = " << exp_limit_val << " cm^2 " << endl;
  cout << "\t +1sigma = " << exp_limit_val + exp_limit_high_1s << endl;
  cout << "\t -1sigma = " << exp_limit_val - exp_limit_low_1s << endl;
  cout << "Observed median limit for mass " << massPoint << "GeV  = " << obs_limit_val << " cm^2 " << endl;
	


//write to file 
   
TFile f("limits/limit_"+massPoint+".root","RECREATE");

obs_limits->Write("obs_limits");
Exp_limits->Write("Exp_limits");
Exp_limitsS1->Write("Exp_limitsS1");
Exp_limitsS2->Write("Exp_limitsS2");

f.Close();
 cout << "END  loop ---------------- " << endl;



/*
///----------------------  HERE START PULLS
  RooAbsPdf *pdftmp = mc.GetPdf();					//get the pdf

  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);			//		
  RooFitResult *fitres = pdftmp->fitTo( data , Save() , Minos(true), PrintLevel(-1) ); //Do the fit on pseudo data

    const RooArgSet *ParaGlobalFit = mc.GetNuisanceParameters();
    RooRealVar* poi = (RooRealVar*) mc.GetParametersOfInterest()->first();
    double muhat = poi->getVal();
   
    cout << "HERE muhat " << muhat << endl; 

    // PLotting the nuisance paramaters pulls 
    RooRealVar* var = NULL;
    const int nPara = ParaGlobalFit->getSize();

    TH1D *h1Dpull = new TH1D("h1Dpull","Pulls",nPara,0.5, nPara + 0.5);
    TH1D *h1sigma = new TH1D("h1sigma","Pulls",nPara, 0.5 , nPara + 0.5);
    h1Dpull->GetYaxis()->SetRangeUser(-4,4);
    h1Dpull->GetYaxis()->SetTitle("(NP_{fit} - NP_{0}) / #DeltaNP");

    double ib=0;
    TIterator* it2 = mc.GetNuisanceParameters()->createIterator();
    while( (var = (RooRealVar*) it2->Next()) ){
      
      // Not consider nuisance parameter being not associated to syst
      TString varname = var->GetName();
      
      double pull  = var->getVal()  ; // GetValue() return value in unit of sigma
      double errorHi = var->getErrorHi() ; 
      double errorLo = var->getErrorLo() ; 
      

      ib++;
      h1Dpull->GetXaxis()->SetBinLabel(ib,varname);
      h1Dpull->SetBinContent(ib,pull);
      h1Dpull->SetBinError(ib,errorHi);

      h1sigma->SetBinContent(ib,0.);	
      h1sigma->SetBinError(ib,1.);	
  
      } 
    h1Dpull->GetXaxis()->LabelsOption("v"); //Ale
    h1Dpull->Draw();
    h1sigma->SetFillColor(3);
    h1sigma->Draw("same E3");
    h1Dpull->Draw("same");

///----  HERE ENDS
*/
}
