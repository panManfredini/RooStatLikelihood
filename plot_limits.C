{

TFile f("limits/limit_all.root");

TGraphErrors *obs_limits = (TGraphErrors*) f.Get("obs_limits"); 
TGraphErrors *Exp_limits = (TGraphErrors*) f.Get("Exp_limits"); 
TGraphAsymmErrors *Exp_limitsS1 = (TGraphAsymmErrors*) f.Get("Exp_limitsS1");
TGraphAsymmErrors *Exp_limitsS2 = (TGraphAsymmErrors*) f.Get("Exp_limitsS2");


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

Exp_limitsS2->GetYaxis()->SetTitle("#sigma");

Exp_limitsS2->GetXaxis()->SetTitle("M  [GeV]");


Exp_limitsS2->GetXaxis()->SetLimits(9.,5000.);
Exp_limitsS2->GetYaxis()->SetRangeUser(1E-38,5E-36);

Exp_limits->GetXaxis()->SetLimits(9.,5000.);
Exp_limits->GetYaxis()->SetRangeUser(1E-38,5E-36);

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

Exp_limits_xmass->SetFillColor(0);
Exp_limits_xmass->SetMarkerSize(0);
Exp_limits_xmass->SetLineWidth(3);
Exp_limits_xmass->SetLineColor(2);


Exp_limitsS2->SetTitle("");

Exp_limitsS2->Draw("AE3");
Exp_limitsS1->Draw("sameE3");
Exp_limits_xmass->Draw("samePc");
Exp_limits->Draw("samePc");
obs_limits->SetMarkerStyle(20);
obs_limits->Draw("samePc");


TLegend* lego = new TLegend(0.2,0.9,0.5,0.7);
  lego->SetTextSize(0.033);
  lego->SetFillColor(0);
  lego->SetBorderSize(0);
  lego->AddEntry(obs_limits,"XENON100 Observed 90\% CLs limit");
  lego->AddEntry(Exp_limits, "XENON100 Expected 90\% CLs limit");
  lego->AddEntry(Exp_limits_xmass, "XMASS 90\% CL limit");
  lego->AddEntry(Exp_limitsS1,"1 #sigma","f");
  lego->AddEntry(Exp_limitsS2,"2 #sigma","f");
  lego->Draw();


gPad->SetLogy();
gPad->SetLogx();
gPad->RedrawAxis("g");

}


