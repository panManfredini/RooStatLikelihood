void limit_reb()
{
//=========Macro generated from canvas: limits/limit
//=========  (Thu Sep 22 12:17:10 2016) by ROOT version6.02/08
   TCanvas *limits = new TCanvas("limits", "limit",550,195,702,600);
   limits->Range(0.6111515,-38.33737,4.042061,-34.96366);
   limits->SetFillColor(0);
   limits->SetBorderMode(0);
   limits->SetBorderSize(2);
   limits->SetLogx();
   limits->SetLogy();
   limits->SetFrameBorderMode(0);
   limits->SetFrameBorderMode(0);
   
   Double_t Graph0_fx3001[18] = {
   20,
   30,
   40,
   50,
   60,
   70,
   80,
   90,
   100,
   200,
   300,
   400,
   500,
   700,
   1000,
   2000,
   3000,
   5000};
   Double_t Graph0_fy3001[18] = {
   6.235417e-36,
   6.656746e-37,
   2.262722e-37,
   1.202159e-37,
   8.13632e-38,
   6.332249e-38,
   5.372319e-38,
   4.833697e-38,
   4.489232e-38,
   4.332769e-38,
   5.252924e-38,
   6.342704e-38,
   7.477645e-38,
   9.802517e-38,
   1.341753e-37,
   2.543024e-37,
   3.764757e-37,
   6.161271e-37};
   Double_t Graph0_felx3001[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fely3001[18] = {
   3.098445e-36,
   3.311771e-37,
   1.125534e-37,
   5.97395e-38,
   4.033772e-38,
   3.144153e-38,
   2.672203e-38,
   2.40628e-38,
   2.232825e-38,
   2.156849e-38,
   2.620711e-38,
   3.168511e-38,
   3.736759e-38,
   4.903225e-38,
   6.71478e-38,
   1.272391e-37,
   1.884008e-37,
   3.081215e-37};
   Double_t Graph0_fehx3001[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fehy3001[18] = {
   6.738838e-36,
   7.334375e-37,
   2.460543e-37,
   1.290625e-37,
   8.677873e-38,
   6.732755e-38,
   5.736414e-38,
   5.168521e-38,
   4.806038e-38,
   4.698138e-38,
   5.711198e-38,
   6.923562e-38,
   8.160256e-38,
   1.072586e-37,
   1.472118e-37,
   2.791597e-37,
   4.134292e-37,
   6.739957e-37};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(18,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
   grae->SetName("Graph0");
   grae->SetTitle("");
   grae->SetFillColor(5);
   grae->SetLineColor(5);
   grae->SetMarkerColor(5);
   grae->SetMarkerSize(0);
   
   TH1F *Graph_Graph3001 = new TH1F("Graph_Graph3001","",100,9,5000);
   Graph_Graph3001->SetMinimum(1e-38);
   Graph_Graph3001->SetMaximum(5e-36);
   Graph_Graph3001->SetDirectory(0);
   Graph_Graph3001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph3001->SetLineColor(ci);
   Graph_Graph3001->GetXaxis()->SetTitle("M_{ #chi}  [GeV]");
   Graph_Graph3001->GetXaxis()->SetLabelFont(42);
   Graph_Graph3001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetXaxis()->SetTitleFont(42);
   Graph_Graph3001->GetYaxis()->SetTitle("SD WIMP Inelastic cross section   [cm^{2}]");
   Graph_Graph3001->GetYaxis()->SetLabelFont(42);
   Graph_Graph3001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetYaxis()->SetTitleOffset(1.3);
   Graph_Graph3001->GetYaxis()->SetTitleFont(42);
   Graph_Graph3001->GetZaxis()->SetLabelFont(42);
   Graph_Graph3001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3001->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3001);
   
   grae->Draw("ae3");
   
   Double_t Graph1_fx3002[18] = {
   20,
   30,
   40,
   50,
   60,
   70,
   80,
   90,
   100,
   200,
   300,
   400,
   500,
   700,
   1000,
   2000,
   3000,
   5000};
   Double_t Graph1_fy3002[18] = {
   6.235417e-36,
   6.656746e-37,
   2.262722e-37,
   1.202159e-37,
   8.13632e-38,
   6.332249e-38,
   5.372319e-38,
   4.833697e-38,
   4.489232e-38,
   4.332769e-38,
   5.252924e-38,
   6.342704e-38,
   7.477645e-38,
   9.802517e-38,
   1.341753e-37,
   2.543024e-37,
   3.764757e-37,
   6.161271e-37};
   Double_t Graph1_felx3002[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fely3002[18] = {
   1.919487e-36,
   2.054087e-37,
   6.969767e-38,
   3.708937e-38,
   2.497175e-38,
   1.945069e-38,
   1.652413e-38,
   1.49063e-38,
   1.382935e-38,
   1.342629e-38,
   1.630138e-38,
   1.96447e-38,
   2.311742e-38,
   3.024679e-38,
   4.142171e-38,
   7.849003e-38,
   1.161671e-37,
   1.90178e-37};
   Double_t Graph1_fehx3002[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fehy3002[18] = {
   2.91616e-36,
   3.144374e-37,
   1.060046e-37,
   5.600252e-38,
   3.780837e-38,
   2.930568e-38,
   2.503852e-38,
   2.248217e-38,
   2.089524e-38,
   2.035206e-38,
   2.469941e-38,
   2.98816e-38,
   3.527862e-38,
   4.641931e-38,
   6.36456e-38,
   1.207218e-37,
   1.788054e-37,
   2.918083e-37};
   grae = new TGraphAsymmErrors(18,Graph1_fx3002,Graph1_fy3002,Graph1_felx3002,Graph1_fehx3002,Graph1_fely3002,Graph1_fehy3002);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(3);
   grae->SetLineColor(3);
   grae->SetMarkerColor(3);
   grae->SetMarkerSize(0);
   
   TH1F *Graph_Graph3002 = new TH1F("Graph_Graph3002","Graph",100,18,5498);
   Graph_Graph3002->SetMinimum(2.691127e-38);
   Graph_Graph3002->SetMaximum(1.006374e-35);
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
   
   grae->Draw("e3");
   
   Double_t Graph2_fx1001[16] = {
   20,
   29.8071,
   39.90202,
   53.41583,
   62.16429,
   69.85718,
   83.21777,
   90,
   105.0887,
   200,
   300,
   388.2045,
   590.8438,
   746.1269,
   1000,
   4244.204};
   Double_t Graph2_fy1001[16] = {
   8e-36,
   7.162923e-37,
   2.027528e-37,
   9.91722e-38,
   7.461589e-38,
   6.3506e-38,
   5.354015e-38,
   5e-38,
   4.600252e-38,
   5e-38,
   6e-38,
   7.252295e-38,
   9.823615e-38,
   1.210266e-37,
   1.5e-37,
   4.354065e-37};
   Double_t Graph2_fex1001[16] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph2_fey1001[16] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphErrors *gre = new TGraphErrors(16,Graph2_fx1001,Graph2_fy1001,Graph2_fex1001,Graph2_fey1001);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");
   gre->SetLineColor(2);
   gre->SetLineWidth(3);
   gre->SetMarkerSize(0);
   
   TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,18,4666.624);
   Graph_Graph1001->SetMinimum(4.140227e-38);
   Graph_Graph1001->SetMaximum(8.7954e-36);
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
   
   gre->Draw("pc");
   
   Double_t Graph3_fx1002[18] = {
   20,
   30,
   40,
   50,
   60,
   70,
   80,
   90,
   100,
   200,
   300,
   400,
   500,
   700,
   1000,
   2000,
   3000,
   5000};
   Double_t Graph3_fy1002[18] = {
   6.235417e-36,
   6.656746e-37,
   2.262722e-37,
   1.202159e-37,
   8.13632e-38,
   6.332249e-38,
   5.372319e-38,
   4.833697e-38,
   4.489232e-38,
   4.332769e-38,
   5.252924e-38,
   6.342704e-38,
   7.477645e-38,
   9.802517e-38,
   1.341753e-37,
   2.543024e-37,
   3.764757e-37,
   6.161271e-37};
   Double_t Graph3_fex1002[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph3_fey1002[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(18,Graph3_fx1002,Graph3_fy1002,Graph3_fex1002,Graph3_fey1002);
   gre->SetName("Graph3");
   gre->SetTitle("Graph");
   gre->SetLineStyle(7);
   gre->SetLineWidth(3);
   gre->SetMarkerSize(0);
   
   TH1F *Graph_Graph1002 = new TH1F("Graph_Graph1002","Graph",100,9,5000);
   Graph_Graph1002->SetMinimum(1e-38);
   Graph_Graph1002->SetMaximum(5e-36);
   Graph_Graph1002->SetDirectory(0);
   Graph_Graph1002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1002->SetLineColor(ci);
   Graph_Graph1002->GetXaxis()->SetLabelFont(42);
   Graph_Graph1002->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1002->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1002->GetXaxis()->SetTitleFont(42);
   Graph_Graph1002->GetYaxis()->SetLabelFont(42);
   Graph_Graph1002->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1002->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1002->GetYaxis()->SetTitleFont(42);
   Graph_Graph1002->GetZaxis()->SetLabelFont(42);
   Graph_Graph1002->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1002->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1002);
   
   gre->Draw("pc");
   
   Double_t Graph4_fx1003[18] = {
   20,
   30,
   40,
   50,
   60,
   70,
   80,
   90,
   100,
   200,
   300,
   400,
   500,
   700,
   1000,
   2000,
   3000,
   5000};
   Double_t Graph4_fy1003[18] = {
   6.799454e-36,
   6.289303e-37,
   1.870549e-37,
   9.204725e-38,
   5.949289e-38,
   4.510389e-38,
   3.771558e-38,
   3.355162e-38,
   3.097464e-38,
   2.942263e-38,
   3.569344e-38,
   4.311977e-38,
   5.08662e-38,
   6.667747e-38,
   9.131451e-38,
   1.729417e-37,
   2.564344e-37,
   4.191275e-37};
   Double_t Graph4_fex1003[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph4_fey1003[18] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(18,Graph4_fx1003,Graph4_fy1003,Graph4_fex1003,Graph4_fey1003);
   gre->SetName("Graph4");
   gre->SetTitle("Graph");
   gre->SetLineWidth(3);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(0);
   
   TH1F *Graph_Graph1003 = new TH1F("Graph_Graph1003","Graph",100,18,5498);
   Graph_Graph1003->SetMinimum(2.648037e-38);
   Graph_Graph1003->SetMaximum(7.476458e-36);
   Graph_Graph1003->SetDirectory(0);
   Graph_Graph1003->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1003->SetLineColor(ci);
   Graph_Graph1003->GetXaxis()->SetLabelFont(42);
   Graph_Graph1003->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph1003->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1003->GetXaxis()->SetTitleFont(42);
   Graph_Graph1003->GetYaxis()->SetLabelFont(42);
   Graph_Graph1003->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph1003->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph1003->GetYaxis()->SetTitleFont(42);
   Graph_Graph1003->GetZaxis()->SetLabelFont(42);
   Graph_Graph1003->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1003->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1003);
   
   gre->Draw("pc");
   
   TLegend *leg = new TLegend(0.307047,0.6660839,0.6073826,0.8653846,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.033);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph4","XENON100 Observed 90% CLs limit","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph3","XENON100 Expected 90% CLs limit","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(7);
   entry->SetLineWidth(3);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","XMASS 90% CL limit","lpf");
   entry->SetFillStyle(1001);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
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
   
   TH1F *Graph_copy = new TH1F("Graph_copy","",100,9,5000);
   Graph_copy->SetMinimum(1e-38);
   Graph_copy->SetMaximum(5e-36);
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
   limits->Modified();
   limits->cd();
   limits->SetSelected(limits);
}
