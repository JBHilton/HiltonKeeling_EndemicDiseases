The scripts in this repository can produce all of the results presented in Hilton and Keeling's paper "Incorporating household structure and demography into models of endemic disease".

ClassedMixingUK and ClassedMixingKenya calculate the age-structured transmission matrices for the UK- and Kenya-like demographic settings, the "translation" matrix E, along with the proportion of each population in each age class at demographic equilibrium.

HH_demo_structured calculates the transition matrices Q_demo, Q_int, and Q_ext, along with their equilibrium H_Eq, the equilibrium prevalence I_bar, the demographic-class stratified prevalence I_T, the class distribution H_T, and childhood infection probability P_R.

Add_demography calculates the demographic transition matrix Q_demo, and updates Q_int to cover all demographic classes.

Grab_Epi_Data calculates early growth parameters R_* and r.

FINISHED EDITS:
ClassedMixingUK
ClassedMixingKenya
HH_demo_structured
Add_Demography
Grab_Epi_Data
State_to_Class
State_to_Class_Mat
State_to_Class_Mat_With_Elders