{
   ifstream in;
   in.open("integral_mass.dat");
   

  vector <double> masses_v;
  vector <double> expected_gaudenz;
  vector <double> expected_xmass;

  double mass_itr =0.;
  double Nev_exp_th_itr =0.;
  double xsec_modifier = 0.;
  double xsec_modifier_xmass = 0.;

  while(mass_itr <1000.){
	in >> mass_itr;
	in >> Nev_exp_th_itr;

 	
	xsec_modifier =  7e-38 *  1.37590955945e-05 / Nev_exp_th_itr; 
	xsec_modifier_xmass = 1.08e-37 * 5.7776674123e-06 / Nev_exp_th_itr;

	masses_v.push_back(mass_itr);
	expected_gaudenz.push_back( xsec_modifier );
	expected_xmass.push_back(xsec_modifier_xmass);
   }


in.close();

const int n = masses_v.size();
double mA[n];
double expected[n];
double expected_v_xmass[n];

for(int k=0; k< n; k++){

	mA[k] = masses_v[k];
	expected[k] = expected_gaudenz[k];
	expected_v_xmass[k] = expected_xmass[k];
}

//double expected_xmass[16] = {8e-36,7e-37, 2e-37, 1e-37, 8e-38, 6e-38, 5.5e-38, 5e-38,  4.3e-38, 5e-38, 6e-38, 7e-38, 9e-38, 1.2e-37, 1.5e-37, 3e-37};
//double m_xmass[16] = { 20, 30., 40., 50., 60., 70., 80., 90.,  100., 200., 300., 400, 500.,700., 1000. ,5000.};


double m_xmass_gg[2] = { 10,5000.};
double m_xmass_g[2] = { 1e-42,1e-33.};

TGraphErrors *Exp_limits_gaudenz = new TGraphErrors(n,mA,expected  );

TGraphErrors *Exp_limits_xmass = new TGraphErrors(n, mA, expected_v_xmass );
TGraphErrors *Exp_limits_xmass_gg = new TGraphErrors(2, m_xmass_gg, m_xmass_g );
Exp_limits_xmass_gg->Draw("AP");
Exp_limits_xmass->Draw("PC");

Exp_limits_gaudenz->Draw("PC");

}
