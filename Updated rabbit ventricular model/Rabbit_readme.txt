Updated version of the Shannon et al. rabbit ventricular model
(Biophys J. 2004 Nov;87(5):3351-71).

This model describes excitation-contraction coupling in the rabbit ventricular myocyte
with integrated descriptions of Ca, Ca/calmodulin-dependent protein kinase II (CaMKII)
and protein kinase A signaling pathways.

The current version was developed by updating the Bartos et al. rabbit ventricular model
(J Physiol. 2017 Apr 1;595(7):2253-2268).

______________________________________________________________________________________
Contents:

Rabbit_readme.txt			        	this file

RabbitVentricularMyocyte_masterCompute.m	loads initial conditions and runs the simulation
RabbitVentricularMyocyte_masterODEfile.m	integrates the following model components
RabbitVentricularMyocyte_barODEfile.m		beta-adrenergic (PKA) phosphorylation module
RabbitVentricularMyocyte_camkiiODEfile.m	CaMKII phosphorylation module
RabbitVentricularMyocyte_camODEfile.m		CaM module
RabbitVentricularMyocyte_eccODEfile.m		excitation-contraction coupling module

.mat files					initial conditions (obtained at 1-Hz pacing)
- yf_rvm_1Hz_NEW				WT model
- yf_rvm_1Hz_120s_ISO_NEW			WT model + 120 s ISO administration (0.1 uM)

- AP_Bartos					voltage signal for AP-clamp simulations
______________________________________________________________________________________


Reference:

S. Morotti, C. Liu, B. Hegyi, H. Ni, A. Fogli Iseppe, L. Wang, M. Pritoni,
C.M. Ripplinger, D.M. Bers, A.G. Edwards, E. Grandi.
Quantitative cross-species translators of cardiac myocyte electrophysiology:
model training, experimental validation, and applications.
Science Advances. 2021 Nov 19;7(47):eabg0927.

Please cite the above paper when using this code.
