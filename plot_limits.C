{


gROOT->ProcessLine(".L ../XeMICRO/Inelastic/inelastic_style.C "); 
setInelasticStyle();

TCanvas canv("canv", "First canvas", 1000, 800);
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.04);

TFile f("limits/smooth_limit_all.root");
//TFile f("limits/limit_all.root");

TGraphErrors *obs_limits = (TGraphErrors*) f.Get("obs_limits"); 
TGraphErrors *Exp_limits = (TGraphErrors*) f.Get("Exp_limits"); 
TGraphAsymmErrors *Exp_limitsS1 = (TGraphAsymmErrors*) f.Get("Exp_limitsS1");
TGraphAsymmErrors *Exp_limitsS2 = (TGraphAsymmErrors*) f.Get("Exp_limitsS2");



Exp_limitsS1->SetFillColor(425);
Exp_limitsS1->SetLineColor(425);
Exp_limitsS1->SetMarkerColor(425);
/*Exp_limitsS1->SetFillColor(851);
Exp_limitsS1->SetLineColor(851);
Exp_limitsS1->SetMarkerColor(851);*/
Exp_limitsS1->SetMarkerSize(0);

Exp_limitsS2->SetFillColor(861);
Exp_limitsS2->SetLineColor(861);
Exp_limitsS2->SetMarkerColor(861);
/*Exp_limitsS2->SetFillColor(852);
Exp_limitsS2->SetLineColor(852);
Exp_limitsS2->SetMarkerColor(852);*/
Exp_limitsS2->SetMarkerSize(0);

obs_limits->SetFillColor(0);
obs_limits->SetLineColor(601);
obs_limits->SetLineWidth(5);
obs_limits->SetMarkerSize(0);

Exp_limits->SetFillColor(0);
Exp_limits->SetMarkerSize(0);
Exp_limits->SetLineStyle(7);
Exp_limits->SetLineWidth(3);

Exp_limitsS2->GetYaxis()->SetTitle("SD Inelastic WIMP cross section  [cm^{2}]");
Exp_limitsS2->GetYaxis()->SetTitleOffset(1.4);

Exp_limitsS2->GetXaxis()->SetTitle("M_{#chi}  [GeV]");


Exp_limitsS2->GetXaxis()->SetLimits(9.,5000.);
Exp_limitsS2->GetYaxis()->SetRangeUser(1E-38,5E-36);

//Exp_limits->GetXaxis()->SetLimits(9.,5000.);
//Exp_limits->GetYaxis()->SetRangeUser(1E-38,5E-36);


double x[16];
double y[16];
x[0]=20; y[0]=8.11629e-36;
x[1]=23; y[1]=2.95962e-36;
x[2]=30; y[2]=6.79689e-37;
x[3]=40; y[3]=2.12121e-37;
x[4]=50; y[4]=1.11727e-37;
x[5]=70; y[5]=6.08576e-38;
x[6]=100; y[6]=4.57758e-38;
x[7]=130; y[7]=4.36316e-38;
x[8]=160; y[8]=4.48429e-38;
x[9]=200; y[9]=4.76474e-38;
x[10]=300; y[10]=6.07221e-38;
x[11]=600; y[11]=1.02375e-37;
x[12]=1000; y[12]=1.60589e-37;
x[13]=2000; y[13]=3.05561e-37;
x[14]=3000; y[14]=4.57505e-37;
x[15]=5000; y[15]=7.46806e-37;

TGraph *Exp_limits_xmass = new TGraph(16, x,y);

Exp_limits_xmass->SetFillColor(0);
Exp_limits_xmass->SetMarkerSize(0);
Exp_limits_xmass->SetMarkerColor(603);
Exp_limits_xmass->SetLineWidth(5);
Exp_limits_xmass->SetLineColor(603);
Exp_limits_xmass->SetLineStyle(9);


Exp_limitsS2->SetTitle("");

//Exp_limits->Draw("APc");
Exp_limitsS2->Draw("A4");
Exp_limitsS1->Draw("same4");
Exp_limits_xmass->Draw("samePc");
obs_limits->SetMarkerStyle(20);
obs_limits->Draw("samePc");


TLegend* lego = new TLegend(0.35,0.90,0.6,0.7);
  lego->SetTextFont(132);
  lego->SetTextSize(0.05);
  lego->SetFillColor(0);
  lego->SetBorderSize(0);
  lego->AddEntry(Exp_limits_xmass, "XMASS   90\% CL limit");
  lego->AddEntry(obs_limits,"XENON100   90\% CLs limit");
//  lego->AddEntry(Exp_limits, "XENON100 Expected 90\% CLs limit");
  lego->AddEntry(Exp_limitsS1,"1 #sigma","f");
  lego->AddEntry(Exp_limitsS2,"2 #sigma","f");
  lego->Draw();


gPad->SetLogy();
gPad->SetLogx();
gPad->RedrawAxis("g");

  canv.Print("/home/pan/work/XeMICRO/Inelastic/SignalSamples/limit_reb.png");
  canv.Print("/home/pan/work/XeMICRO/Inelastic/SignalSamples/limit_reb.eps");
  canv.Print("/home/pan/Documents/Papers/mypapers/Inelastic/limit_reb.png");
  canv.Print("/home/pan/Documents/Papers/mypapers/Inelastic/limit_reb.eps");


}


