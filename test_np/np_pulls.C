#define np_pulls_cxx
#include "np_pulls.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TProfile.h>

void np_pulls::Loop()
{
//   In a ROOT session, you can do:
//      root> .L np_pulls.C
//      root> np_pulls t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

  vector <double*> nuissance;
nuissance.push_back(&T_Acc);
nuissance.push_back(&T_Bkg_norm_factor);
nuissance.push_back(&T_eff_1);
nuissance.push_back(&T_eff_10);
nuissance.push_back(&T_eff_11);
nuissance.push_back(&T_eff_12);
nuissance.push_back(&T_eff_13);
nuissance.push_back(&T_eff_14);
nuissance.push_back(&T_eff_15);
nuissance.push_back(&T_eff_16);
nuissance.push_back(&T_eff_17);
nuissance.push_back(&T_eff_18);
nuissance.push_back(&T_eff_19);
nuissance.push_back(&T_eff_2);
nuissance.push_back(&T_eff_20);
nuissance.push_back(&T_eff_3);
nuissance.push_back(&T_eff_4);
nuissance.push_back(&T_eff_5);
nuissance.push_back(&T_eff_6);
nuissance.push_back(&T_eff_7);
nuissance.push_back(&T_eff_8);
nuissance.push_back(&T_eff_9);
nuissance.push_back(&T_n_Co_1);
nuissance.push_back(&T_n_Co_10);
nuissance.push_back(&T_n_Co_11);
nuissance.push_back(&T_n_Co_12);
nuissance.push_back(&T_n_Co_13);
nuissance.push_back(&T_n_Co_14);
nuissance.push_back(&T_n_Co_15);
nuissance.push_back(&T_n_Co_16);
nuissance.push_back(&T_n_Co_17);
nuissance.push_back(&T_n_Co_18);
nuissance.push_back(&T_n_Co_19);
nuissance.push_back(&T_n_Co_2);
nuissance.push_back(&T_n_Co_20);
nuissance.push_back(&T_n_Co_3);
nuissance.push_back(&T_n_Co_4);
nuissance.push_back(&T_n_Co_5);
nuissance.push_back(&T_n_Co_6);
nuissance.push_back(&T_n_Co_7);
nuissance.push_back(&T_n_Co_8);
nuissance.push_back(&T_n_Co_9);
nuissance.push_back(&mu);

   TH2F *correlation_matrix = new TH2F("correlation_matrix","correlation matrix", 50,0.5,50.5,50,0.5,50.5);
   correlation_matrix->GetXaxis()->SetBinLabel(1,"T_Acc");
   correlation_matrix->GetXaxis()->SetBinLabel(2,"T_Bkg_norm_factor");
   correlation_matrix->GetXaxis()->SetBinLabel(3, "T_eff_1");
   correlation_matrix->GetXaxis()->SetBinLabel(4,"T_eff_10");
   correlation_matrix->GetXaxis()->SetBinLabel(5,"T_eff_11");
   correlation_matrix->GetXaxis()->SetBinLabel(6,"T_eff_12");
   correlation_matrix->GetXaxis()->SetBinLabel(7,"T_eff_13");
   correlation_matrix->GetXaxis()->SetBinLabel(8,"T_eff_14");
   correlation_matrix->GetXaxis()->SetBinLabel(9,"T_eff_15");
   correlation_matrix->GetXaxis()->SetBinLabel(10,"T_eff_16");
   correlation_matrix->GetXaxis()->SetBinLabel(11,"T_eff_17");
   correlation_matrix->GetXaxis()->SetBinLabel(12,"T_eff_18");
   correlation_matrix->GetXaxis()->SetBinLabel(13,"T_eff_19");
   correlation_matrix->GetXaxis()->SetBinLabel(14,"T_eff_2");
   correlation_matrix->GetXaxis()->SetBinLabel(15,"T_eff_20");
   correlation_matrix->GetXaxis()->SetBinLabel(16,"T_eff_3");
   correlation_matrix->GetXaxis()->SetBinLabel(17,"T_eff_4");
   correlation_matrix->GetXaxis()->SetBinLabel(18,"T_eff_5");
   correlation_matrix->GetXaxis()->SetBinLabel(19,"T_eff_6");
   correlation_matrix->GetXaxis()->SetBinLabel(20,"T_eff_7");
   correlation_matrix->GetXaxis()->SetBinLabel(21,"T_eff_8");
   correlation_matrix->GetXaxis()->SetBinLabel(22,"T_eff_9");
   correlation_matrix->GetXaxis()->SetBinLabel(23,"T_n_Co_1");
   correlation_matrix->GetXaxis()->SetBinLabel(24,"T_n_Co_10");
   correlation_matrix->GetXaxis()->SetBinLabel(25,"T_n_Co_11");
   correlation_matrix->GetXaxis()->SetBinLabel(26,"T_n_Co_12");
   correlation_matrix->GetXaxis()->SetBinLabel(27,"T_n_Co_13");
   correlation_matrix->GetXaxis()->SetBinLabel(28,"T_n_Co_14");
   correlation_matrix->GetXaxis()->SetBinLabel(29,"T_n_Co_15");
   correlation_matrix->GetXaxis()->SetBinLabel(30,"T_n_Co_16");
   correlation_matrix->GetXaxis()->SetBinLabel(31,"T_n_Co_17");
   correlation_matrix->GetXaxis()->SetBinLabel(32,"T_n_Co_18");
   correlation_matrix->GetXaxis()->SetBinLabel(33,"T_n_Co_19");
   correlation_matrix->GetXaxis()->SetBinLabel(34,"T_n_Co_2");
   correlation_matrix->GetXaxis()->SetBinLabel(35,"T_n_Co_20");
   correlation_matrix->GetXaxis()->SetBinLabel(36,"T_n_Co_3");
   correlation_matrix->GetXaxis()->SetBinLabel(37,"T_n_Co_4");
   correlation_matrix->GetXaxis()->SetBinLabel(38,"T_n_Co_5");
   correlation_matrix->GetXaxis()->SetBinLabel(39,"T_n_Co_6");
   correlation_matrix->GetXaxis()->SetBinLabel(40,"T_n_Co_7");
   correlation_matrix->GetXaxis()->SetBinLabel(41,"T_n_Co_8");
   correlation_matrix->GetXaxis()->SetBinLabel(42,"T_n_Co_9");
   correlation_matrix->GetXaxis()->SetBinLabel(43,"mu");

   correlation_matrix->GetYaxis()->SetBinLabel(1,"T_Acc");
   correlation_matrix->GetYaxis()->SetBinLabel(2,"T_Bkg_norm_factor");
   correlation_matrix->GetYaxis()->SetBinLabel(3, "T_eff_1");
   correlation_matrix->GetYaxis()->SetBinLabel(4,"T_eff_10");
   correlation_matrix->GetYaxis()->SetBinLabel(5,"T_eff_11");
   correlation_matrix->GetYaxis()->SetBinLabel(6,"T_eff_12");
   correlation_matrix->GetYaxis()->SetBinLabel(7,"T_eff_13");
   correlation_matrix->GetYaxis()->SetBinLabel(8,"T_eff_14");
   correlation_matrix->GetYaxis()->SetBinLabel(9,"T_eff_15");
   correlation_matrix->GetYaxis()->SetBinLabel(10,"T_eff_16");
   correlation_matrix->GetYaxis()->SetBinLabel(11,"T_eff_17");
   correlation_matrix->GetYaxis()->SetBinLabel(12,"T_eff_18");
   correlation_matrix->GetYaxis()->SetBinLabel(13,"T_eff_19");
   correlation_matrix->GetYaxis()->SetBinLabel(14,"T_eff_2");
   correlation_matrix->GetYaxis()->SetBinLabel(15,"T_eff_20");
   correlation_matrix->GetYaxis()->SetBinLabel(16,"T_eff_3");
   correlation_matrix->GetYaxis()->SetBinLabel(17,"T_eff_4");
   correlation_matrix->GetYaxis()->SetBinLabel(18,"T_eff_5");
   correlation_matrix->GetYaxis()->SetBinLabel(19,"T_eff_6");
   correlation_matrix->GetYaxis()->SetBinLabel(20,"T_eff_7");
   correlation_matrix->GetYaxis()->SetBinLabel(21,"T_eff_8");
   correlation_matrix->GetYaxis()->SetBinLabel(22,"T_eff_9");
   correlation_matrix->GetYaxis()->SetBinLabel(23,"T_n_Co_1");
   correlation_matrix->GetYaxis()->SetBinLabel(24,"T_n_Co_10");
   correlation_matrix->GetYaxis()->SetBinLabel(25,"T_n_Co_11");
   correlation_matrix->GetYaxis()->SetBinLabel(26,"T_n_Co_12");
   correlation_matrix->GetYaxis()->SetBinLabel(27,"T_n_Co_13");
   correlation_matrix->GetYaxis()->SetBinLabel(28,"T_n_Co_14");
   correlation_matrix->GetYaxis()->SetBinLabel(29,"T_n_Co_15");
   correlation_matrix->GetYaxis()->SetBinLabel(30,"T_n_Co_16");
   correlation_matrix->GetYaxis()->SetBinLabel(31,"T_n_Co_17");
   correlation_matrix->GetYaxis()->SetBinLabel(32,"T_n_Co_18");
   correlation_matrix->GetYaxis()->SetBinLabel(33,"T_n_Co_19");
   correlation_matrix->GetYaxis()->SetBinLabel(34,"T_n_Co_2");
   correlation_matrix->GetYaxis()->SetBinLabel(35,"T_n_Co_20");
   correlation_matrix->GetYaxis()->SetBinLabel(36,"T_n_Co_3");
   correlation_matrix->GetYaxis()->SetBinLabel(37,"T_n_Co_4");
   correlation_matrix->GetYaxis()->SetBinLabel(38,"T_n_Co_5");
   correlation_matrix->GetYaxis()->SetBinLabel(39,"T_n_Co_6");
   correlation_matrix->GetYaxis()->SetBinLabel(40,"T_n_Co_7");
   correlation_matrix->GetYaxis()->SetBinLabel(41,"T_n_Co_8");
   correlation_matrix->GetYaxis()->SetBinLabel(42,"T_n_Co_9");
   correlation_matrix->GetYaxis()->SetBinLabel(43,"mu");

   TProfile *nuis_pulls  = new TProfile("nuis_pulls","Pulls",50,0.5,50.5,-10.,10.,"s");
   nuis_pulls->GetXaxis()->SetBinLabel(1,"T_Acc");
   nuis_pulls->GetXaxis()->SetBinLabel(2,"T_Bkg_norm_factor");
   nuis_pulls->GetXaxis()->SetBinLabel(3, "T_eff_1");
   nuis_pulls->GetXaxis()->SetBinLabel(4,"T_eff_10");
   nuis_pulls->GetXaxis()->SetBinLabel(5,"T_eff_11");
   nuis_pulls->GetXaxis()->SetBinLabel(6,"T_eff_12");
   nuis_pulls->GetXaxis()->SetBinLabel(7,"T_eff_13");
   nuis_pulls->GetXaxis()->SetBinLabel(8,"T_eff_14");
   nuis_pulls->GetXaxis()->SetBinLabel(9,"T_eff_15");
   nuis_pulls->GetXaxis()->SetBinLabel(10,"T_eff_16");
   nuis_pulls->GetXaxis()->SetBinLabel(11,"T_eff_17");
   nuis_pulls->GetXaxis()->SetBinLabel(12,"T_eff_18");
   nuis_pulls->GetXaxis()->SetBinLabel(13,"T_eff_19");
   nuis_pulls->GetXaxis()->SetBinLabel(14,"T_eff_2");
   nuis_pulls->GetXaxis()->SetBinLabel(15,"T_eff_20");
   nuis_pulls->GetXaxis()->SetBinLabel(16,"T_eff_3");
   nuis_pulls->GetXaxis()->SetBinLabel(17,"T_eff_4");
   nuis_pulls->GetXaxis()->SetBinLabel(18,"T_eff_5");
   nuis_pulls->GetXaxis()->SetBinLabel(19,"T_eff_6");
   nuis_pulls->GetXaxis()->SetBinLabel(20,"T_eff_7");
   nuis_pulls->GetXaxis()->SetBinLabel(21,"T_eff_8");
   nuis_pulls->GetXaxis()->SetBinLabel(22,"T_eff_9");
   nuis_pulls->GetXaxis()->SetBinLabel(23,"T_n_Co_1");
   nuis_pulls->GetXaxis()->SetBinLabel(24,"T_n_Co_10");
   nuis_pulls->GetXaxis()->SetBinLabel(25,"T_n_Co_11");
   nuis_pulls->GetXaxis()->SetBinLabel(26,"T_n_Co_12");
   nuis_pulls->GetXaxis()->SetBinLabel(27,"T_n_Co_13");
   nuis_pulls->GetXaxis()->SetBinLabel(28,"T_n_Co_14");
   nuis_pulls->GetXaxis()->SetBinLabel(29,"T_n_Co_15");
   nuis_pulls->GetXaxis()->SetBinLabel(30,"T_n_Co_16");
   nuis_pulls->GetXaxis()->SetBinLabel(31,"T_n_Co_17");
   nuis_pulls->GetXaxis()->SetBinLabel(32,"T_n_Co_18");
   nuis_pulls->GetXaxis()->SetBinLabel(33,"T_n_Co_19");
   nuis_pulls->GetXaxis()->SetBinLabel(34,"T_n_Co_2");
   nuis_pulls->GetXaxis()->SetBinLabel(35,"T_n_Co_20");
   nuis_pulls->GetXaxis()->SetBinLabel(36,"T_n_Co_3");
   nuis_pulls->GetXaxis()->SetBinLabel(37,"T_n_Co_4");
   nuis_pulls->GetXaxis()->SetBinLabel(38,"T_n_Co_5");
   nuis_pulls->GetXaxis()->SetBinLabel(39,"T_n_Co_6");
   nuis_pulls->GetXaxis()->SetBinLabel(40,"T_n_Co_7");
   nuis_pulls->GetXaxis()->SetBinLabel(41,"T_n_Co_8");
   nuis_pulls->GetXaxis()->SetBinLabel(42,"T_n_Co_9");
   nuis_pulls->GetXaxis()->SetBinLabel(43,"mu");

   double xy [50][50] = {0}; //sum of product xy
   double x [50] = {0};	// sum X	
   double variance_np [50] = {0};	// variance of X	
   double x_2 [50] = {0};	// sum X^2

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

  	cout << *nuissance[0] << endl;
 	for(int n=0; n < nuissance.size() ; n++)
	{
	   nuis_pulls->Fill((double)n+1.,*nuissance[n], 1.);
 	   x[n] += *nuissance[n];
	   x_2[n] += pow(*nuissance[n],2.);
 	   
	   for(int k =0; k < nuissance.size() ; k++)
	   {
		xy[n][k] += (*nuissance[n]) * (*nuissance[k]);
	   }  
	}	

   }



   for(int i=0;i <50; i++){

	variance_np[i] = sqrt(x_2[i] / nentries - pow(x[i]/nentries, 2.));  
//	nuis_pulls->SetBinError(i+1, variance_np[i]);
 	cout <<"var " <<  variance_np[i] << endl;

      for(int j=0;j <50; j++){
	double correlation = (nentries*xy[j][i] - x[j]*x[i] ) /sqrt( ( nentries*x_2[j] - pow(x[j],2.)) * (nentries*x_2[i] - pow(x[i],2.)) ) ;
   	correlation_matrix->SetBinContent(j+1,i+1,correlation);
      }
   }


   nuis_pulls->GetYaxis()->SetRangeUser(-1.5,1.5);
   nuis_pulls->Draw();

   correlation_matrix->GetZaxis()->SetRangeUser(-1.,1.);
   new TCanvas();
   correlation_matrix->Draw("colz");


   TFile *f=new TFile(nameFile,"RECREATE");
   f->cd();
   correlation_matrix->Write();
   nuis_pulls->Write();
   f->Close();


}







