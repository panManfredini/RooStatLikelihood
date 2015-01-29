
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
#include "RooRandom.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooMath.h"
#include "Math/QuantFuncMathCore.h"


using namespace RooStats ;
using namespace RooFit ;

void SimplestIntervalExample() {
// set RooFit random seed for reproducible results
RooRandom::randomGenerator()->SetSeed(3001) ;
// make a simple model via the workspace factory
RooWorkspace * wspace = new RooWorkspace () ;
wspace->factory("Gaussian::normal(x[-10,10],mu[ -1 ,1],sigma[1])") ;
wspace->defineSet( "poi" ,"mu" ) ;
wspace->defineSet( "obs" ,"x" ) ;
// specify components of model for statistical tools
ModelConfig *modelConfig = new ModelConfig("Example G(x|mu,1)" ) ;
modelConfig-> SetWorkspace (*wspace ) ;
modelConfig->SetPdf( *wspace->pdf( "normal" )) ;
modelConfig->SetParametersOfInterest(*wspace->set("poi")) ;
modelConfig->SetObservables(*wspace->set( "obs")) ;
// create a toy dataset
RooDataSet * data = wspace->pdf( "normal" )->generate(*wspace->set("obs") ,10) ;
data->Print() ;
// for convenience later on
RooRealVar *x = wspace-> var( "x" ) ;
RooRealVar *mu = wspace->var( "mu" ) ;
// set confidence level


double conf = 0.95;
// example use profile likelihood calculator
ProfileLikelihoodCalculator plc (*data , *modelConfig ) ;
plc.SetConfidenceLevel(conf) ;
LikelihoodInterval *plInt = plc.GetInterval() ;
// some plots
TCanvas * canvas = new TCanvas ( "canvas" , "canvas" , 600 ,300) ;
canvas->Divide(2) ;
// plot the data
canvas->cd(1) ;
RooPlot * frame = x->frame() ;
data->plotOn(frame) ;
data->statOn(frame) ;
frame->Draw() ;

// plot the profile likeihood
canvas->cd(2) ;
LikelihoodIntervalPlot plot(plInt) ;
plot.Draw()  ;
// for this example we know the expected intervals analytically
double expectedLL = data->mean(*x) + ROOT::Math::normal_quantile( (1 - conf)/2 ,1)/ sqrt ( data->numEntries() ) ;
double expectedUL = data->mean(*x) + ROOT::Math::normal_quantile_c((1 -conf) /2 ,1)/ sqrt ( data->numEntries() ) ;
// Print the
std:: cout << endl << " expected interval is [ " << expectedLL <<  " , " << expectedUL <<  " ] " << endl ;

cout << " plc interval is [ " << plInt->LowerLimit(*mu) << " , " << plInt->UpperLimit (*mu) << " ] " << endl ;

}

