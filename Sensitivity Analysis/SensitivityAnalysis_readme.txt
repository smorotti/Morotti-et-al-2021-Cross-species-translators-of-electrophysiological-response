Codes required to perform parameter sensitivity analysis on the human, rabbit,
and mouse models paced at 1 Hz in control condition (Fig. 1 and Fig. S2)

______________________________________________________________________________________
PART 1 - Generation of the simulated population-based dataset
(note that the commands for saving the results are now commented)

1) Run SA_ALL_01_generate_parameters to generate the matrix with the scaling factors
	(1000 model variants created with sigma = 0.1)

2) Run the following codes to generate and analyze the population of human models:
	SA_ALL_Human_02_obtain_ICs
	SA_ALL_Human_03_beat_analysis

3) Run the following codes to generate and analyze the population of rabbit models:
	SA_ALL_Rabbit_02_obtain_ICs
	SA_ALL_Rabbit_03_beat_analysis
    
4) Run the following codes to generate and analyze the population of mouse models:
	SA_ALL_Mouse_02_obtain_ICs
	SA_ALL_Mouse_03_beat_analysis

______________________________________________________________________________________
PART 2 - Sensitivity analysis
    
5) Run "SA_ALL_04_comparison_PLS_analysis" to perform analysis and plot results

______________________________________________________________________________________


Reference:

S. Morotti, C. Liu, B. Hegyi, H. Ni, A. Fogli Iseppe, L. Wang, M. Pritoni,
C.M. Ripplinger, D.M. Bers, A.G. Edwards, E. Grandi.
Quantitative cross-species translators of cardiac myocyte electrophysiology:
model training, experimental validation, and applications.
Science Advances. 2021 Nov 19;7(47):eabg0927.

Please cite the above paper when using these codes.
