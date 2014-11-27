#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"

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


using namespace RooStats ;
using namespace RooFit ;

void test_counting_experiment() {

//////////////////////   MODEL BUILDING    /////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/*
int nobs = 3;   // number of observed events
double b0 = 1; // number of background events
double sigmab = 0.2;   // relative uncertainty in b
*/

RooWorkspace w("w");

w.factory("prod:s(eff[1,0,50], Xsec[100.], K_m[20], mu[1,0,200])");
w.factory("sum:nexp(s, b[50,0,200])");
// Poisson of (n | s+b)
w.factory("Poisson:pdf(nobs[0,200],nexp)");
// Poisson constraint to express uncertainty in b
w.factory("Poisson:constraint_B(b0[50,0,200],b)");  // Likelihood for nuissance parameter b
w.factory("prod:N_S(S_tot[50], eff)");
w.factory("Poisson:constraint_S(S0[50,0,200], N_S)");  // Likelihood for nuissance parameter b
// the total model p.d.f will be the product of the two
//w.factory("Uniform::prior_b(b)"); // define the prior for b
//w.factory("Uniform::prior_s(s)"); // define the prior for b

w.factory("PROD:model(pdf,constraint_S,constraint_B)");


w.var("b0")->setConstant(true); 
w.var("S_tot")->setConstant(true); 
w.var("S0")->setConstant(true); 
w.var("Xsec")->setConstant(true); 
w.var("K_m")->setConstant(true); 

ModelConfig mc("ModelConfig",&w);
mc.SetPdf(*w.pdf("model"));
mc.SetParametersOfInterest(*w.var("mu"));
mc.SetObservables(*w.var("nobs"));

//mc.SetPriorPdf(*w.pdf("prior_s"));
mc.SetNuisanceParameters(*w.var("b,eff"));

// need now to set the global observable
mc.SetGlobalObservables(*w.var("b0"));

// this is needed for the hypothesis tests
mc.SetSnapshot(*w.var("mu"));

// import model in the workspace 
w.import(mc);

// make data set with the number of observed events
RooDataSet data("data","", *w.var("nobs"));
w.var("nobs")->setVal(55);
data.add(*w.var("nobs") );

// import data set in workspace and save it in a file
w.import(data);
w.writeToFile("CountingModel.root", true);

w.Print();

data.Print();


/*///////////////////////////////////////////////////////////////////////
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
  plot->SetRange(0,50);  // possible eventually to change ranges
  //plot->SetNPoints(50);  // do not use too many points, it could become very slow for some models
  plot->Draw("");  // use option TF1 if too slow (plot.Draw("tf1")


*/
///-------------------- bayesian calcularot ----//
/*
  BayesianCalculator bayesianCalc(data,mc);
  bayesianCalc.SetConfidenceLevel(0.95); // 68% interval

  // set the type of interval (not really needed for central which is the default)
  //bayesianCalc.SetLeftSideTailFraction(0.5); // for central interval
  bayesianCalc.SetLeftSideTailFraction(0.); // for upper limit
  //bayesianCalc.SetShortestInterval(); // for shortest interval


  // set the integration type (not really needed for the default ADAPTIVE)
  // possible alternative values are  "VEGAS" , "MISER", or "PLAIN"  (MC integration from libMathMore) 
  // "TOYMC" (toy MC integration, work when nuisances exist and they have a constraints pdf)
  TString integrationType = "TOYMC";

  // this is needed if using TOYMC
  if (integrationType.Contains("TOYMC") ) { 
    RooAbsPdf * nuisPdf = RooStats::MakeNuisancePdf(mc, "constraint");
    if (nuisPdf) bayesianCalc.ForceNuisancePdf(*nuisPdf);
  }

  // limit the number of iterations to speed-up the integral 
  bayesianCalc.SetNumIters(10000);
  bayesianCalc.SetIntegrationType(integrationType); 

  // compute interval by scanning the posterior function
  // it is done by default when computing shortest intervals
  bayesianCalc.SetScanOfPosterior(50);

  RooRealVar* firstPOI = (RooRealVar*) mc.GetParametersOfInterest()->first();

  SimpleInterval* interval = bayesianCalc.GetInterval();
  if (!interval) { 
     cout << "Error computing Bayesian interval - exit " << endl;
     return;
  }


  double lowerLimit = interval->LowerLimit();
  double upperLimit = interval->UpperLimit();


  cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
    lowerLimit << ", "<<
    upperLimit <<"] "<<endl;

  // draw plot of posterior function

  RooPlot * plot = bayesianCalc.GetPosteriorPlot();
  if (plot) plot->Draw();  
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


//------------ Getting the interval as function of m --------------//
  double massMin = 10.;
  double massMax = 100.;

  vector <double> masses_v;
  vector <double> observed_v;
  vector <double> expected_v;
  vector <double> expected_S1_up_v;
  vector <double> expected_S1_dw_v;
  vector <double> expected_S2_up_v;
  vector <double> expected_S2_dw_v;

  for( double mass=massMin; mass<=massMax; mass += (massMax-massMin)/10.0 ){

	  w.var("K_m")->setVal(1./ mass);

	  AsymptoticCalculator  ac(data, *bModel, *sbModel);
	  //ac.SetOneSidedDiscovery(true);  // for one-side discovery test
	  ac.SetOneSided(true);  // for one-side tests (limits)
	  //  ac->SetQTilde(true);
	  //ac.SetPrintLevel(-1);  // to suppress print level 


	// create hypotest inverter 
	  // passing the desired calculator 
	  HypoTestInverter *calc = new HypoTestInverter(ac);    // for asymptotic 
	  //HypoTestInverter calc(fc);  // for frequentist

	  calc->SetConfidenceLevel(0.95);
	  calc->UseCLs(true);
	  int npoints = 500;  // number of points to scan
	  // min and max (better to choose smaller intervals)
	  double poimin = poi->getMin();
	  double poimax = poi->getMax();
	  //poimin = 0; poimax=10;

	  std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
	  calc->SetFixedScan(npoints,poimin,poimax);
  
	  HypoTestInverterResult * r = calc->GetInterval();

	  double upperLimit = r->UpperLimit();

/*	  std::cout << "The computed upper limit is: " << upperLimit << std::endl;
  	  // compute expected limit
  	std::cout << "Expected upper limits, using the B (alternate) model on mu: " << std::endl;
  	std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
  	std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
  	std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
  	std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
  	std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;
  */
  	std::cout << "Expected upper limits, using the B (alternate) model in terms of Xsec: " << std::endl;
  	std::cout << "The computed upper limit is: " << w.var("Xsec")->getValV() * upperLimit << std::endl;
  	std::cout << " expected limit (median) " << w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(0) << std::endl;
  	std::cout << " expected limit (-1 sig) " << w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(-1) << std::endl;
  	std::cout << " expected limit (+1 sig) " << w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(1) << std::endl;
  	std::cout << " expected limit (-2 sig) " << w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(-2) << std::endl;
  	std::cout << " expected limit (+2 sig) " << w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(2) << std::endl;
	
	masses_v.push_back(mass);
	observed_v.push_back( w.var("Xsec")->getValV() * upperLimit );
	expected_v.push_back( w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(0) );
	expected_S1_up_v.push_back(w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(1));
	expected_S2_up_v.push_back(w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(2));
	expected_S2_dw_v.push_back(w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(-2));
	expected_S1_dw_v.push_back(w.var("Xsec")->getValV() * r->GetExpectedUpperLimit(-1));
	
/*	observed_v.push_back( w.var("Xsec")->getValV() *  w.var("K_m")->getValV()* upperLimit );
	expected_v.push_back( w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(0) );
	expected_S1_up_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(1));
	expected_S2_up_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(2));
	expected_S2_dw_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(-2));
	expected_S1_dw_v.push_back(w.var("Xsec")->getValV() * w.var("K_m")->getValV()* r->GetExpectedUpperLimit(-1));
*/
   }


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
//obs_limits->SetMarkerStyle(5);

Exp_limits->SetFillColor(0);
Exp_limits->SetMarkerSize(0);
Exp_limits->SetLineStyle(7);
Exp_limits->SetLineWidth(3);


//Exp_limitsS2->GetYaxis()->SetTitle("#sigma#timesBR( #phi #rightarrow #tau#tau )  [pb]");
Exp_limitsS2->GetYaxis()->SetTitle("#sigma");

Exp_limitsS2->GetXaxis()->SetTitle("M  [GeV]");

Exp_limitsS2->Draw("Al3");
Exp_limitsS1->Draw("sameL3");
Exp_limits->Draw("PL");
obs_limits->Draw("PL");

Exp_limitsS2->GetYaxis()->SetRangeUser(0.,100.);

TLegend* lego = new TLegend(0.2,0.9,0.5,0.7);
  lego->SetTextSize(0.033);
  lego->SetFillColor(0);
  lego->SetBorderSize(0);
  lego->AddEntry(obs_limits,"Observed 95\% CLs limit");
  lego->AddEntry(Exp_limits, "Expected 95\% CLs limit");
  lego->AddEntry(Exp_limitsS1,"1 #sigma","f");
  lego->AddEntry(Exp_limitsS2,"2 #sigma","f");
  lego->Draw();


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











