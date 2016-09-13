{

 ifstream in;
 in.open("pppp");

 double test =0.;
 double integral =0;

 for(int i=0; i<200; i++){
	in >> test;
	integral += test;

	cout << test << endl;
 }
  
  cout << "Integral " << integral << endl;
  cout << "Rate with abundance " << integral * 0.2644 << endl;

}
