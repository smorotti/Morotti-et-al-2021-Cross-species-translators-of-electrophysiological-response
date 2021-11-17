Updated version of the Morotti et al. mouse ventricular model
(J Physiol. 2014 Mar 15;592(6):1181-97).

This model describes excitation-contraction coupling in the mouse ventricular myocyte
with integrated descriptions of Ca, Ca/calmodulin-dependent protein kinase II (CaMKII)
and protein kinase A signaling pathways.

The current version was developed by updating the Surdo et al. mouse ventricular model
(Nat Commun. 2017 Apr 20;8:15031).

______________________________________________________________________________________
Contents:

Mouse_readme.txt			        this file

MouseVentricularMyocyte_masterCompute.m	loads initial conditions and runs the simulation
MouseVentricularMyocyte_masterODEfile.m	integrates the following model components
MouseVentricularMyocyte_barODEfile.m	beta-adrenergic (PKA) phosphorylation module
MouseVentricularMyocyte_camkiiODEfile.m	CaMKII phosphorylation module
MouseVentricularMyocyte_camODEfile.m	CaM module
MouseVentricularMyocyte_eccODEfile.m	excitation-contraction coupling module

.mat files				initial conditions (obtained at 1-Hz pacing)
- yf_mvm_1Hz_NEW			WT model
- yf_mvm_1Hz_120s_ISO_NEW		WT model + 120 s ISO administration (0.1 uM)
______________________________________________________________________________________


Reference:

S. Morotti, C. Liu, B. Hegyi, H. Ni, A. Fogli Iseppe, L. Wang, M. Pritoni,
C.M. Ripplinger, D.M. Bers, A.G. Edwards, E. Grandi.
Quantitative cross-species translators of cardiac myocyte electrophysiology:
model training, experimental validation, and applications.
Science Advances. 2021 Nov 19;7(47):eabg0927.

Please cite the above paper when using this code.
