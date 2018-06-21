{


TFile f("limits/limit_all.root");

TGraphErrors *obs_limits = (TGraphErrors*) f.Get("obs_limits"); 
TGraphErrors *Exp_limits = (TGraphErrors*) f.Get("Exp_limits"); 
TGraphAsymmErrors *Exp_limitsS1 = (TGraphAsymmErrors*) f.Get("Exp_limitsS1");
TGraphAsymmErrors *Exp_limitsS2 = (TGraphAsymmErrors*) f.Get("Exp_limitsS2");

Double_t* Ey_high = Exp_limitsS1->GetEYhigh();
Double_t* Ey_low = Exp_limitsS1->GetEYlow();

Double_t* Ey_high_2s = Exp_limitsS2->GetEYhigh();
Double_t* Ey_low_2s = Exp_limitsS2->GetEYlow();



Double_t* x = Exp_limitsS1->GetX();
Double_t* y = Exp_limitsS1->GetY();

TGraph median;
TGraph Es1_high;
TGraph Es2_high;
TGraph Es1_low;
TGraph Es2_low;

for(int i=0; i< Exp_limitsS1->GetN(); i++){

	median.SetPoint(i, x[i], y[i]);
	Es1_high.SetPoint(i, x[i], Ey_high[i] + y[i]);	
	Es2_high.SetPoint(i, x[i], Ey_high_2s[i] + y[i]);	
	Es1_low.SetPoint(i, x[i], y[i] - Ey_low[i]);
	Es2_low.SetPoint(i, x[i], y[i] - Ey_low_2s[i]);

	cout << i << "   "  << x[i] << endl;
}


median.Draw("APC");
Es1_high.Draw("samePC");
Es2_high.Draw("samePC");
Es1_low.Draw("samePC");
Es2_low.Draw("samePC");


TGraphAsymmErrors *smoothed_Exp_limitsS1 = new TGraphAsymmErrors();
TGraphAsymmErrors *smoothed_Exp_limitsS2 = new TGraphAsymmErrors();


for(i =0 ; i <=5 ; i++){

	smoothed_Exp_limitsS1->SetPoint(i, x[i], y[i]);
	smoothed_Exp_limitsS1->SetPointError(i, 0., 0., Ey_low[i], Ey_high[i]  );
	smoothed_Exp_limitsS2->SetPoint(i, x[i], y[i]);
	smoothed_Exp_limitsS2->SetPointError(i, 0., 0.,  Ey_low_2s[i], Ey_high_2s[i]  );

}

for(j =6 ; j < 10000; j++){
	
	double X = 80. + ((double) j) * 5000. / 10000.;
        double Y = median.Eval(X, 0, "S");
	double es1_low = Y - Es1_low.Eval(X, 0, "S") ; 	
	double es2_low = Y - Es2_low.Eval(X, 0, "S") ; 	
	double es1_high = Es1_high.Eval(X, 0, "S") - Y;
	double es2_high = Es2_high.Eval(X, 0, "S") - Y;

	smoothed_Exp_limitsS1->SetPoint(j, X, Y);
	smoothed_Exp_limitsS1->SetPointError(j, 0., 0., es1_low, es1_high );
	smoothed_Exp_limitsS2->SetPoint(j, X, Y);
	smoothed_Exp_limitsS2->SetPointError(j, 0., 0., es2_low, es2_high );

}


smoothed_Exp_limitsS1->SetFillColor(4);
smoothed_Exp_limitsS2->SetFillColor(3);
smoothed_Exp_limitsS2->Draw("same3");
smoothed_Exp_limitsS1->Draw("same3");



TFile *foo = new TFile("limits/smooth_limit_all.root","RECREATE");

foo->cd();


obs_limits->Write("obs_limits");
Exp_limits->Write("Exp_limits");
smoothed_Exp_limitsS1->Write("Exp_limitsS1");
smoothed_Exp_limitsS2->Write("Exp_limitsS2");


foo->Close();

}












