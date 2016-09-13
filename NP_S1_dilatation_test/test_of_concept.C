{

gStyle->SetOptStat(0);

// test per np che mimica una S1 = Gauss * multiplier
// che succede se prendo una gauss e la moltiplico? come posso codificare questo effetto su un histo?

TH1F *h = new TH1F("h","",50,0.,10.);
TH1F *h2 = new TH1F("h2","",50,0.,10.);

double a = 1.05; // multiplicative factor

TF1 gaus1 ("gaus1","[0]*exp(-0.5*((x/[3]-[1])/[2])**2)/[3]",0.,10.);
gaus1.SetParameters(1.,5.,1.,1.);

h->FillRandom("gaus1",1000000);

TCanvas *c1 = new TCanvas();
h->Draw();
gaus1.SetParameter(3, a);
h2->FillRandom("gaus1",1000000);
h2->SetLineColor(3);
h2->Draw("same");

TCanvas *c2 = new TCanvas();

TH1F *h_clone = (TH1F*) h2->Clone("ciccio");

h_clone->Divide(h);

h_clone->Draw();

TF1 factor ("factor","exp(0.5*(x^2*(1- 1/([2]**2))/([1]^2)))* exp(x*[0]*(-1. + 1./[2])/([1]**2))/[2]",0.,10.);
factor.SetParameters(5.,1.,a);

factor.Draw("same");


TH1F *h_clone2 = (TH1F*) h->Clone("ciccio2");

for(int j=1; j <=50; j++){
	double x = h_clone2->GetBinCenter(j);
	h_clone2->SetBinContent(j, h_clone2->GetBinContent(j) * factor.Eval(x));
   }

c1->cd();
//new TCanvas();
h_clone2->SetLineColor(6);
//h_clone2->Divide(h2);
//h_clone2->Draw();
h_clone2->Draw("same");
}
