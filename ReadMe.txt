The scripts in this repository can produce all of the results presented in Hilton and Keeling's paper "Incorporating household structure and demography into models of endemic disease".

______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
NOTES ON THE SCRIPTS AND FUNCTIONS:

State_to_Class calculates the expected number of age class C individuals in a household which is in demographic class T. State_to_Class_Mat calculates this for all age classes when there is no "elder" class - State_to_Class_Mat_With_Elders does the same but with an "elder" class.

ClassedMixingUK and ClassedMixingKenya calculate the age-structured transmission matrices for the UK- and Kenya-like demographic settings, the "translation" matrix E, along with the proportion of each population in each age class at demographic equilibrium.

HH_demo_structured calculates the transition matrices Q_demo, Q_int, and Q_ext, along with their equilibrium H_Eq, the equilibrium prevalence I_bar, the demographic-class stratified prevalence I_T, the class distribution H_T, and childhood infection probability P_R.

Add_demography calculates the demographic transition matrix Q_demo, and updates Q_int to cover all demographic classes.

Grab_Epi_Data calculates early growth parameters R_* and r.

GetHouseholdSizeDistributions calculates the disease-free equilibria of the two demographic settings and plots their household size distributions, separated into young and old households. These plots are in Figure 2 of the paper.

The two Get[]Demography scripts calculate the demographic parameters used in all the loop scripts so that we don't need to repeatedly calculate them.

UKStructure comparison calculates the behaviour of mumps in the UK under the four transmission structures. This is the data plotted in Figure 3 of the paper and Figure 1 of the supplement, and listed in Table 2 of the paper.

Get_Prevalence_Histogram is a function which takes the equilibrium distribution and plots the prevalence stratified by household demography, as seen in Figure 3 of the paper etc.

Get_StructureHists plots the histograms from Figure 3 of the paper and Figure 1 of the supplement.

The four SigmaLoop[][] scripts perform sensitivity analysis by calculating early growth parameters and equilibrium behaviour of the model over varying values of the tuning parameter sigma. Data generated by this script are plotted in Figure 4 of the paper and Figures 2 and 9-12 of the supplement, and listed in Table 3 of the paper.

GetPinkBarData gets the data for the homogeneous baseline comparison plotted as pink bars in Figure 4 of the paper and Figure 2 of the supplement.

Get_DiseaseDemogHists plots the histograms from Figure 4 of the paper and Figure 2 of the supplement.

PlotSigmaResults plots the sensitivity analysis results from Figures 9-12 of the supplement.

PlotAgeBurdens plots the histograms from Figures 3 and 4 of the supplement.

The four TauLoop[][] scripts perform sensitivity analysis by calculating early growth parameters and equilibrium behaviour of the model over varying values of the transmission rate tau. Data generated by this script are plotted in Figures 5-8 of the supplement.

PlotTauResults plots the sensitivity analysis results from Figures 5-8 of the supplement.

PlotContactMatrices plots the heatmaps From Figure 14 of the supplement.

______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
NOTES ON THE DATA:

ContactData.zip contains contact survey data used to calculate the contact durations used in the model. Parameters.zip contains all the parameters we use in the calculations presented in the paper. ModelOuput.zip contains all the outputs presented in the paper. When the scripts in the repository are run, outputs are saved to one of these folders, depending on the nature of the output.

All of the .mat files calculated by these scripts have a time and date label in their filename. The 'original' outputs calculated by the authors are provided in Parameters.zip and ModelOutput.zip and have no such label in their filename.

The two .mat files in ContactData.zip contain data derived from the supplementary information attached to the papers by Kiti et. al. and Mossong et. al.  cited in our paper.