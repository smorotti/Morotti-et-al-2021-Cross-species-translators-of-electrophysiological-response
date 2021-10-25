Updated version of the Grandi et al. human ventricular model
(J Mol Cell Cardiol. 2010 Jan;48(1):112-21).

This model describes excitation-contraction coupling in the human ventricular myocyte
with integrated descriptions of Ca, Ca/calmodulin-dependent protein kinase II (CaMKII)
and protein kinase A signaling pathways.

The current version was developed by updating the Moreno et al. human ventricular model
(Circ Res. 2013 Sep 13;113(7):e50-e61).

______________________________________________________________________________________
Contents:

Human_readme.txt			        	this file

HumanVentricularMyocyte_masterCompute.m		loads initial conditions and runs the simulation
HumanVentricularMyocyte_masterODEfile.m		integrates the following model components
HumanVentricularMyocyte_barODEfile.m		beta-adrenergic (PKA) phosphorylation module
HumanVentricularMyocyte_camkiiODEfile.m		CaMKII phosphorylation module
HumanVentricularMyocyte_camODEfile.m		CaM module
HumanVentricularMyocyte_eccODEfile.m		excitation-contraction coupling module

.mat files					initial conditions (obtained at 1-Hz pacing)
- yf_hvm_1Hz_NEW				WT model
- yf_hvm_1Hz_120s_ISO_NEW			WT model + 120 s ISO administration (0.1 uM)
______________________________________________________________________________________


Reference:

S. Morotti, C. Liu, B. Hegyi, H. Ni, A. Fogli Iseppe, L. Wang, M. Pritoni,
C.M. Ripplinger, D.M. Bers, A.G. Edwards, E. Grandi.
Quantitative cross-species translators of cardiac myocyte electrophysiology:
model training, experimental validation, and applications.
In revision for Sci Adv.

Please cite the above paper when using this code.
