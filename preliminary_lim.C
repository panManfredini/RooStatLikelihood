void preliminary_lim()
{
//=========Macro generated from canvas: limits/limit
//=========  (Thu Apr  7 17:02:10 2016) by ROOT version6.02/08
   TCanvas *limits = new TCanvas("limits", "limit",447,166,600,600);
   limits->Range(0.6985228,-38.625,3.25572,-32.375);
   limits->SetFillColor(0);
   limits->SetBorderMode(0);
   limits->SetBorderSize(2);
   limits->SetLogx();
   limits->SetLogy();
   limits->SetFrameBorderMode(0);
   limits->SetFrameBorderMode(0);
   
   Double_t Graph0_fx3001[8] = {
   20.09233,
   33.51603,
   55.9081,
   93.26033,
   155.5676,
   247.7076,
   413.2012,
   689.2612};
   Double_t Graph0_fy3001[8] = {
   4.564978e-36,
   3.382083e-37,
   8.201299e-38,
   4.325023e-38,
   3.89107e-38,
   4.578761e-38,
   6.373184e-38,
   9.712285e-38};
   Double_t Graph0_felx3001[8] = {
   6.952695e-310,
   6.952688e-310,
   6.953148e-310,
   6.953148e-310,
   6.953148e-310,
   6.952979e-310,
   6.952979e-310,
   6.952979e-310};
   Double_t Graph0_fely3001[8] = {
   2.279653e-36,
   1.682728e-37,
   4.085733e-38,
   2.15569e-38,
   1.938976e-38,
   2.281441e-38,
   3.175766e-38,
   4.839642e-38};
   Double_t Graph0_fehx3001[8] = {
   6.952695e-310,
   6.952688e-310,
   6.953148e-310,
   6.953148e-310,
   6.953148e-310,
   6.952979e-310,
   6.952979e-310,
   6.952979e-310};
   Double_t Graph0_fehy3001[8] = {
   4.95633e-36,
   3.630511e-37,
   8.821099e-38,
   4.683932e-38,
   4.218429e-38,
   4.962391e-38,
   6.904579e-38,
   1.05213e-37};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(8,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillColor(5);
   grae->SetLineColor(5);
   grae->SetMarkerColor(5);
   grae->SetMarkerSize(0);
   
   TH1F *Graph_Graph3001 = new TH1F("Graph_Graph3001","Graph",100,9,1000);
   Graph_Graph3001->SetMinimum(1e-38);
   Graph_Graph3001->SetMaximum(1e-33);
   Graph_Graph3001->SetDirectory(0);
   Graph_Graph3001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph3001->SetLineColor(ci);
   Graph_Graph3001->GetXaxis()->SetTitle("M  [GeV]");
   Graph_Graph3001->GetXaxis()->SetLabelFont(42);
   Graph_Graph3001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetXaxis()->SetTitleFont(42);
   Graph_Graph3001->GetYaxis()->SetTitle("#sigma");
   Graph_Graph3001->GetYaxis()->SetLabelFont(42);
   Graph_Graph3001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetYaxis()->SetTitleFont(42);
   Graph_Graph3001->GetZaxis()->SetLabelFont(42);
   Graph_Graph3001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3001);
   
   grae->Draw("al3");
   
   Double_t Graph1_fx3002[8] = {
   20.09233,
   33.51603,
   55.9081,
   93.26033,
   155.5676,
   247.7076,
   413.2012,
   689.2612};
   Double_t Graph1_fy3002[8] = {
   4.564978e-36,
   3.382083e-37,
   8.201299e-38,
   4.325023e-38,
   3.89107e-38,
   4.578761e-38,
   6.373184e-38,
   9.712285e-38};
   Double_t Graph1_felx3002[8] = {
   6.952695e-310,
   6.952688e-310,
   6.953148e-310,
   6.953148e-310,
   6.953148e-310,
   6.952979e-310,
   6.952979e-310,
   6.952979e-310};
   Double_t Graph1_fely3002[8] = {
   1.413276e-36,
   1.045712e-37,
   2.542257e-38,
   1.338911e-38,
   1.203478e-38,
   1.416058e-38,
   1.971677e-38,
   3.004842e-38};
   Double_t Graph1_fehx3002[8] = {
   6.952695e-310,
   6.952688e-310,
   6.953148e-310,
   6.953148e-310,
   6.953148e-310,
   6.952979e-310,
   6.952979e-310,
   6.952979e-310};
   Double_t Graph1_fehy3002[8] = {
   2.149442e-36,
   1.57988e-37,
   3.83695e-38,
   2.034025e-38,
   1.830985e-38,
   2.15401e-38,
   2.99751e-38,
   4.567748e-38};
   grae = new TGraphAsymmErrors(8,Graph1_fx3002,Graph1_fy3002,Graph1_felx3002,Graph1_fehx3002,Graph1_fely3002,Graph1_fehy3002);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(3);
   grae->SetLineColor(3);
   grae->SetMarkerColor(3);
   grae->SetMarkerSize(0);
   
   TH1F *Graph_Graph3002 = new TH1F("Graph_Graph3002","Graph",100,18.0831,756.1781);
   Graph_Graph3002->SetMinimum(2.418832e-38);
   Graph_Graph3002->SetMaximum(7.383174e-36);
   Graph_Graph3002->SetDirectory(0);
   Graph_Graph3002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3002->SetLineColor(ci);
   Graph_Graph3002->GetXaxis()->SetLabelFont(42);
   Graph_Graph3002->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3002->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3002->GetXaxis()->SetTitleFont(42);
   Graph_Graph3002->GetYaxis()->SetLabelFont(42);
   Graph_Graph3002->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3002->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3002->GetYaxis()->SetTitleFont(42);
   Graph_Graph3002->GetZaxis()->SetLabelFont(42);
   Graph_Graph3002->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3002->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3002->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3002);
   
   grae->Draw("l3");
   
   Double_t Graph2_fx1001[8] = {
   20.09233,
   33.51603,
   55.9081,
   93.26033,
   155.5676,
   247.7076,
   413.2012,
   689.2612};
   Double_t Graph2_fy1001[8] = {
   4.564978e-36,
   3.382083e-37,
   8.201299e-38,
   4.325023e-38,
   3.89107e-38,
   4.578761e-38,
   6.373184e-38,
   9.712285e-38};
   Double_t Graph2_fex1001[8] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph2_fey1001[8] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphErrors *gre = new TGraphErrors(8,Graph2_fx1001,Graph2_fy1001,Graph2_fex1001,Graph2_fey1001);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");
   gre->SetLineStyle(7);
   gre->SetLineWidth(3);
   gre->SetMarkerSize(0);
   
   TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,9,1000);
   Graph_Graph1001->SetMinimum(1e-38);
   Graph_Graph1001->SetMaximum(1e-33);
   Graph_Graph1001->SetDirectory(0);
   Graph_Graph1001->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1001->SetLineColor(ci);
   Graph_Graph1001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1001);
   
   gre->Draw("pl");
   
   TLegend *leg = new TLegend(0.3842282,0.6363636,0.6845638,0.8356643,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.033);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","Observed 90% CLs limit","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","Expected 90% CLs limit","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(7);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","1 #sigma","f");
   entry->SetFillColor(3);
   entry->SetFillStyle(1001);
   entry->SetLineColor(3);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph0","2 #sigma","f");
   entry->SetFillColor(5);
   entry->SetFillStyle(1001);
   entry->SetLineColor(5);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   
   TH1F *Graph_copy = new TH1F("Graph_copy","Graph",100,9,1000);
   Graph_copy->SetMinimum(1e-38);
   Graph_copy->SetMaximum(1e-33);
   Graph_copy->SetDirectory(0);
   Graph_copy->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_copy->SetLineColor(ci);
   Graph_copy->GetXaxis()->SetTitle("M  [GeV]");
   Graph_copy->GetXaxis()->SetLabelFont(42);
   Graph_copy->GetXaxis()->SetLabelSize(0.035);
   Graph_copy->GetXaxis()->SetTitleSize(0.035);
   Graph_copy->GetXaxis()->SetTitleFont(42);
   Graph_copy->GetYaxis()->SetTitle("#sigma");
   Graph_copy->GetYaxis()->SetLabelFont(42);
   Graph_copy->GetYaxis()->SetLabelSize(0.035);
   Graph_copy->GetYaxis()->SetTitleSize(0.035);
   Graph_copy->GetYaxis()->SetTitleFont(42);
   Graph_copy->GetZaxis()->SetLabelFont(42);
   Graph_copy->GetZaxis()->SetLabelSize(0.035);
   Graph_copy->GetZaxis()->SetTitleSize(0.035);
   Graph_copy->GetZaxis()->SetTitleFont(42);
   Graph_copy->Draw("sameaxis");
   
   TPaveText *pt = new TPaveText(0.4195302,0.9341608,0.5804698,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *AText = pt->AddText("Graph");
   pt->Draw();
   limits->Modified();
   limits->cd();
   limits->SetSelected(limits);
}
