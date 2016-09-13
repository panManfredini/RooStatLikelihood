{

   ifstream in;
   in.open("integral_mass.dat");

  double mass_itr =0.;
  double Nev_exp_th_itr =0.;
  double xsec_modifier = 10.;
//  double N_tot_theory = w.var("N_tot_theory")->getValV();

  while(mass_itr <1000.){
	in >> mass_itr;
	in >> Nev_exp_th_itr;
	int ciccio = mass_itr;
        cout << ciccio << endl;
  }


}
