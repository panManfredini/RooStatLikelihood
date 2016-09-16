#!bash

echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("20");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_20  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("30");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_30  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("40");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_40  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("50");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_50  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("60");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_60  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("70");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_70  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("80");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_80  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("90");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_90  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("100");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_100  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("200");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_200  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("300");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_300  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("400");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_400  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("500");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_500  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("700");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_700  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("1000");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_1000  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("2000");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_2000  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("3000");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_4000  &
echo 'gROOT->LoadMacro("limit_inelastic_stable.C"); limit_inelastic_stable("5000");  gSystem->Exit(0);'  | root -b -l &> LOGS/log_5000  &

