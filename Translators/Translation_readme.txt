Codes required to build, analyze, and apply the translators of electrophyiological
response across species or experimental condition.

The theoretical methodology applied here was originally developed by Gong and Sobie
(NPJ Syst Biol Appl. 2018 Feb 24;4:11).


Reference:

S. Morotti, C. Liu, B. Hegyi, H. Ni, A. Fogli Iseppe, L. Wang, M. Pritoni,
C.M. Ripplinger, D.M. Bers, A.G. Edwards, E. Grandi.
Quantitative cross-species translators of cardiac myocyte electrophysiology:
model training, experimental validation, and applications.
Science Advances. 2021 Nov 19;7(47):eabg0927.

Please cite the above paper when using these codes.

______________________________________________________________________________________
PART 1 - Generation of the simulated population-based dataset
(note that the commands for saving the results are now commented)

1) Run SA_COMMON_01_generate_parameters to generate the matrix with scaling factors
affecting only parameters common to the mouse, rabbit, and human models
	(1500 model variants created with sigma = 0.1)

2) Run the following codes to generate and analyze the population of human models
(simulated at the indicated pacing rate/condition):
	SA_COMMON_Human_02_obtain_ICs
	SA_COMMON_Human_03_beat_analysis	(1-Hz, control, steady-state)

3) Run the following codes to generate and analyze the population of rabbit models
(simulated at the indicated pacing rate/condition):
	SA_COMMON_Rabbit_02_obtain_ICs
	SA_COMMON_Rabbit_03_beat_analysis	(1-Hz, control, steady-state)

	SA_COMMON_Rabbit_04_obtain_ICs
	SA_COMMON_Rabbit_05_beat_analysis	(0.5-Hz, control, steady-state)

	SA_COMMON_Rabbit_06_obtain_ICs
	SA_COMMON_Rabbit_07_beat_analysis	(2-Hz, control, steady-state)

	SA_COMMON_Rabbit_08_obtain_ICs
	SA_COMMON_Rabbit_09_beat_analysis	(3-Hz, control, steady-state)

	SA_COMMON_Rabbit_10_obtain_ICs_ISO
	SA_COMMON_Rabbit_11_beat_analysis_ISO	(3.5-Hz, ISO,
		60-s simulation starting from steady-state control condition at 2-Hz)
    
4) Run the following codes to generate and analyze the population of mouse models
(simulated at the indicated pacing rate/condition):
	SA_COMMON_Mouse_02_obtain_ICs
	SA_COMMON_Mouse_03_beat_analysis	(1-Hz, control, steady-state)

	SA_COMMON_Mouse_04_obtain_ICs
	SA_COMMON_Mouse_05_beat_analysis	(4.8-Hz, control, steady-state)

	SA_COMMON_Mouse_06_obtain_ICs
	SA_COMMON_Mouse_07_beat_analysis	(6.2-Hz, ISO,
		60-s simulation starting from steady-state control condition at 4.8-Hz)
    
5) Run the following codes to analyze the response of the baseline mouse, rabbit, and
human models to selective ion channel block or ISO administration (1-Hz):
	baseline_models_analysis_mouse
	baseline_models_analysis_rabbit
	baseline_models_analysis_human

______________________________________________________________________________________
PART 2 - Applications to simulated data & analysis

6) Run "translator_MtoH_or_RtoH_x" to plot the results shown in Fig. 3A-B
(change "input_species" for selecting Mouse or Rabbit as input species).

7) Run "translator_MtoH_or_RtoH_varying_nf" to plot the results shown in Fig. 3C.

8) Run "translate_simulated_drug_effect" to plot the results shown in Fig. 4, Fig. S3,
and Fig. S4 (change current_index_SA to select the desired condition, i.e. control or block).

9) Run "recursive_feature_elimination_analysis" to plot the results shown in Fig. 6
(change input_species_index for selecting Mouse or Rabbit as input species).

______________________________________________________________________________________
PART 3 - Applications to experimental data

10) Run "translate_drug_mouse_to_human" to plot results in Fig. 5A.

11) Run "translate_drug_rabbit_to_human" to plot results in Fig. 5A-C and Fig. S7
(change block_index to select the desired current to block).

12) Run "translate_perturbation_across_frequency_rabbit"  to plot results shown in Fig. 7,
Fig. S8, and Fig. 8A (change drug_index for selecting different ion channel blockers or ISO).

13) Run "translate_ISO_mouse_to_rabbit" to plot results shown in Fig. 8B-D
(set n_features_mouse to 2 or 4 to change the number of mouse features considered).

14) Run "translate_SNS_mouse_to_rabbit" to plot results shown in Fig. 9.

______________________________________________________________________________________
