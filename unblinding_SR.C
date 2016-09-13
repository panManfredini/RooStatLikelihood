{


  TFile *histos  = TFile::Open("histo_limits_gaudenz_final.root");

  double Bkg_norm_factor = 0.0212;//0.0378;

  TH1F  *h_co60_SR = (TH1F*)histos->Get("Co-hist");
  TH1F  *h_DM_SR = (TH1F*)histos->Get("DM-hist");
  TH1F  *h_ambe_SR = (TH1F*)histos->Get("hist11");

  h_co60_SR->Scale(Bkg_norm_factor);
  h_co60_SR->Print();
  h_DM_SR->Print();

  cout << "sigma " << (h_DM_SR->Integral() - h_co60_SR->Integral())/sqrt(h_DM_SR->Integral()) << endl;
  h_co60_SR->SetFillColor(2);

  //new TCanvas();

  //new TCanvas();
  h_ambe_SR->Scale(50./h_ambe_SR->Integral());
  h_ambe_SR->Add(h_co60_SR);
  h_ambe_SR->Draw("hist");
  h_co60_SR->Draw("hist same");
  h_DM_SR->Draw("sameE"); 

TLegend* lego = new TLegend(0.2,0.9,0.5,0.7);
  lego->SetTextSize(0.033);
  lego->SetFillColor(0);
  lego->SetBorderSize(0);
  lego->AddEntry(h_co60_SR, "ER Background","f");
  lego->AddEntry(h_ambe_SR, "5000 GeV WIMP (50 events)","f");
  lego->AddEntry(h_DM_SR, "Data","p");
  lego->Draw();
}
