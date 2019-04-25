The scripts in this repository can produce all of the results presented in Hilton and Keeling's paper "Incorporating household structure and demography into models of endemic disease".

State_to_Class calculates the expected number of age class C individuals in a household which is in demographic class T. State_to_Class_Mat calculates this for all age classes when there is no "elder" class - State_to_Class_Mat_With_Elders does the same but with an "elder" class.

ClassedMixingUK and ClassedMixingKenya calculate the age-structured transmission matrices for the UK- and Kenya-like demographic settings, the "translation" matrix E, along with the proportion of each population in each age class at demographic equilibrium.

HH_demo_structured calculates the transition matrices Q_demo, Q_int, and Q_ext, along with their equilibrium H_Eq, the equilibrium prevalence I_bar, the demographic-class stratified prevalence I_T, the class distribution H_T, and childhood infection probability P_R.

Add_demography calculates the demographic transition matrix Q_demo, and updates Q_int to cover all demographic classes.

Grab_Epi_Data calculates early growth parameters R_* and r.

GetHouseholdSizeDistributions calculates the disease-free equilibria of the two demographic settings and plots their household size distributions, separated into young and old households. These plots are in Figure 2 of the paper.

The two Get[]Demography scripts calculate the demographic parameters used in all the loop scripts so that we don't need to repeatedly calculate them.

UKStructure comparison calculates the behaviour of mumps in the UK under the four transmission structures. This is the data plotted in Figure 3 of the paper and Figure 1 of the supplement, and listed in Table 2 of the paper.


FINISHED EDITS:
ClassedMixingUK
ClassedMixingKenya
HH_demo_structured
Add_Demography
Grab_Epi_Data
State_to_Class
State_to_Class_Mat
State_to_Class_Mat_With_Elders
GetHouseholdSizeDistributions
GetUKDemography
GetKenyaDemography
UKStructureComparison
Get_Prevalence_Histogram
Get_StructureHists

CURRENTLY ON SIGMALOOPUKMUMPS