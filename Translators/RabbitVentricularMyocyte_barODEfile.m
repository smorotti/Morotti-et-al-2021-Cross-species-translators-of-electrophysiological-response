function ydot = RabbitVentricularMyocyte_barODEfile(t,y,pin)
% This module describes the beta-adrenergic signaling pathway.

%% State variables

LR=y(1);
LRG=y(2);
RG=y(3);
b1AR_S464=y(4);
b1AR_S301=y(5);
GsaGTPtot=y(6);
GsaGDP=y(7);
Gsby=y(8);
AC_GsaGTP=y(9);
%PDEp=y(10);
PDE3p=y(10);
PDE4p=y(11);
cAMPtot=y(12);
RC_I=y(13);
RCcAMP_I=y(14);
RCcAMPcAMP_I=y(15);
RcAMPcAMP_I=y(16);
PKACI=y(17);
PKACI_PKI=y(18);
RC_II=y(19);
RCcAMP_II=y(20);
RCcAMPcAMP_II=y(21);
RcAMPcAMP_II=y(22);
PKACII=y(23);
PKACII_PKI=y(24);
I1p_PP1=y(25); % output CaMKII
I1ptot=y(26);
PLBp=y(27); % output
PLMp=y(28); % output
LCCap=y(29); % output
LCCbp=y(30); % output
RyRp=y(31); % output
TnIp=y(32); % output
Myop=y(33); % output
KSp=y(34); % output
KRp=y(35); % output
CFTRp=y(36); % output
ClCap=y(37); % output
Top=y(38); % output
K1p=y(39); % output
Nap=y(40); % output

ydot = zeros(size(y));
%% Drug Concentrations

ISO = pin(1); % (uM) isoproterenol concentration - Ltot
FSK = 0; % (uM) forskolin concentration
IBMX = 0; % (uM) IBMX concentration
%% b-AR module

%b1ARtot = 0.00528; % (uM) total b1-AR protein % MOUSE
%b1ARtot = 0.0132; % RAT
b1ARtot = 0.028; % RABBIT

kf_LR           = 1;              % (1/[uM ms]) forward rate for ISO binding to b1AR
kr_LR           = 0.285;          % (1/ms) reverse rate for ISO binding to b1AR
kf_LRG          = 1;              % (1/[uM ms]) forward rate for ISO:b1AR association with Gs
kr_LRG          = 0.062;          % (1/ms) reverse rate for ISO:b1AR association with Gs
kf_RG           = 1;              % (1/[uM ms]) forward rate for b1AR association with Gs
kr_RG           = 33;             % (1/ms) reverse rate for b1AR association with Gs

Gstot           = 3.83;           % (uM) total Gs protein
k_G_act         = 16e-3;          % (1/ms) rate constant for Gs activation
k_G_hyd         = 0.8e-3;         % (1/ms) rate constant for G-protein hydrolysis
k_G_reassoc     = 1.21;           % (1/[uM ms]) rate constant for G-protein reassociation

kf_bARK         = 1.1e-6;         % (1/[uM ms]) forward rate for b1AR phosphorylation by b1ARK
kr_bARK         = 2.2e-6;         % (1/ms) reverse rate for b1AR phosphorylation by b1ARK
kf_PKA          = 3.6e-6;         % (1/[uM ms]) forward rate for b1AR phosphorylation by PKA
kr_PKA          = 2.2e-6;         % (1/ms) reverse rate for b1AR phosphorylation by PKA

b1ARact = b1ARtot - b1AR_S464 - b1AR_S301;
b1AR = b1ARact - LR - LRG - RG;
Gs = Gstot - LRG - RG - Gsby;

dLR = kf_LR*ISO*b1AR - kr_LR*LR + kr_LRG*LRG - kf_LRG*LR*Gs;
dLRG = kf_LRG*LR*Gs - kr_LRG*LRG - k_G_act*LRG;
dRG = kf_RG*b1AR*Gs - kr_RG*RG - k_G_act*RG;

bARK_desens = kf_bARK*(LR+LRG);
bARK_resens = kr_bARK*b1AR_S464;
PKA_desens = kf_PKA*PKACI*b1ARact;
PKA_resens = kr_PKA*b1AR_S301;
db1AR_S464 = bARK_desens - bARK_resens; % ydot(5)
db1AR_S301 = PKA_desens - PKA_resens; % ydot(6)

G_act = k_G_act*(RG+LRG);
G_hyd = k_G_hyd*GsaGTPtot;
G_reassoc = k_G_reassoc*GsaGDP*Gsby;
dGsaGTPtot = G_act - G_hyd; % ydot(7)
dGsaGDP = G_hyd - G_reassoc; % ydot(8)
dGsby = G_act - G_reassoc; % ydot(9)
% end b-AR module
%% cAMP module

%ACtot = 70.57e-3;       % (uM) total adenylyl cyclase % MOUSE
ACtot = 47e-3; % RABBIT

ATP             = 5e3;            % (uM) total ATP
k_AC_basal      = 0.2e-3;         % (1/ms) basal cAMP generation rate by AC
Km_AC_basal     = 1.03e3;         % (uM) basal AC affinity for ATP

Kd_AC_Gsa       = 0.4;            % (uM) Kd for AC association with Gsa
kf_AC_Gsa       = 1;              % (1/[uM ms]) forward rate for AC association with Gsa
kr_AC_Gsa       = Kd_AC_Gsa;      % (1/ms) reverse rate for AC association with Gsa

k_AC_Gsa        = 8.5e-3;         % (1/ms) basal cAMP generation rate by AC:Gsa
Km_AC_Gsa       = 315.0;          % (uM) AC:Gsa affinity for ATP

Kd_AC_FSK       = 44.0;           % (uM) Kd for FSK binding to AC
k_AC_FSK        = 7.3e-3;         % (1/ms) basal cAMP generation rate by AC:FSK
Km_AC_FSK       = 860.0;          % (uM) AC:FSK affinity for ATP

% MOUSE
%PDEtot          = 22.85e-3;       % (uM) total phosphodiesterase
%k_cAMP_PDE      = 5e-3;           % (1/ms) cAMP hydrolysis rate by PDE
%k_cAMP_PDEp     = 2*k_cAMP_PDE;   % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
%Km_PDE_cAMP     = 1.3;            % (uM) PDE affinity for cAMP

% RABBIT
%PDE3tot = 0.036;       % (uM) total phosphodiesterase 
%PDE4tot = 0.036;       % (uM) total phosphodiesterase
    PDE3tot = 0.036*0.75; % PROVA
    PDE4tot = 0.036*0.75; % PROVA
k_cAMP_PDE3 = 3.5e-3;               % k_pde3        [1/ms]
k_cAMP_PDE3p = 2*k_cAMP_PDE3;   % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
Km_PDE3_cAMP = 0.15;             % Km_pde3       [uM]
k_cAMP_PDE4 = 5.0e-3;               % k_pde4        [1/ms]
k_cAMP_PDE4p = 2*k_cAMP_PDE4;   % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
Km_PDE4_cAMP = 1.3;              % Km_pde4       [uM]

Kd_PDE_IBMX     = 30.0;           % (uM) Kd_R2cAMP_C for IBMX binding to PDE
k_PKA_PDE       = 7.5e-3;         % (1/ms) rate constant for PDE phosphorylation by type 1 PKA
k_PP_PDE        = 1.5e-3;         % (1/ms) rate constant for PDE dephosphorylation by phosphatases

cAMP = cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II);
AC = ACtot - AC_GsaGTP;
GsaGTP = GsaGTPtot - AC_GsaGTP;
dAC_GsaGTP = kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP;

AC_FSK = FSK*AC/Kd_AC_FSK;
AC_ACT_BASAL = k_AC_basal*AC*ATP/(Km_AC_basal+ATP);
AC_ACT_GSA = k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP);
AC_ACT_FSK = k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP);

% MOUSE
% %PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX;
% PDE_IBMX = PDEtot*IBMX/(Kd_PDE_IBMX+IBMX);
% PDE = PDEtot - PDE_IBMX - PDEp;
% dPDEp = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp;
% PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP);
% 
% dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT; % ydot(11)

% RABBIT                                % Add constrain on total IBMX?
%PDE3_IBMX = PDE3tot*IBMX/Kd_PDE_IBMX;
PDE3_IBMX = PDE3tot*IBMX/(Kd_PDE_IBMX+IBMX);
PDE3 = PDE3tot - PDE3_IBMX - PDE3p;
dPDE3p = k_PKA_PDE*PKACII*PDE3 - k_PP_PDE*PDE3p; % ydot(10)
PDE3_ACT = k_cAMP_PDE3*PDE3*cAMP/(Km_PDE3_cAMP+cAMP) + k_cAMP_PDE3p*PDE3p*cAMP/(Km_PDE3_cAMP+cAMP);

%PDE4_IBMX = PDE4tot*IBMX/Kd_PDE_IBMX;
PDE4_IBMX = PDE4tot*IBMX/(Kd_PDE_IBMX+IBMX);
PDE4 = PDE4tot - PDE4_IBMX - PDE4p;
dPDE4p = k_PKA_PDE*PKACII*PDE4 - k_PP_PDE*PDE4p;  % ydot() - NEW STATE VARIABLE
PDE4_ACT = k_cAMP_PDE4*PDE4*cAMP/(Km_PDE4_cAMP+cAMP) + k_cAMP_PDE4p*PDE4p*cAMP/(Km_PDE4_cAMP+cAMP);

dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE3_ACT - PDE4_ACT; % ydot(11)
% end cAMP module
%% PKA module

PKItot          = 0.18;           % (uM) total PKI
kf_RC_cAMP      = 1;              % (1/[uM ms]) Kd for PKA RC binding to cAMP
kf_RCcAMP_cAMP  = 1;              % (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
kf_RcAMPcAMP_C  = 4.375;          % (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
kf_PKA_PKI      = 1;              % (1/[uM ms]) Ki for PKA inhibition by PKI
kr_RC_cAMP      = 1.64;           % (1/ms) Kd for PKA RC binding to cAMP
kr_RCcAMP_cAMP  = 9.14;           % (1/ms) Kd for PKA RC:cAMP binding to cAMP
kr_RcAMPcAMP_C  = 1;              % (1/ms) Kd for PKA R:cAMPcAMP binding to C
kr_PKA_PKI      = 2e-4;           % (1/ms) Ki for PKA inhibition by PKI

epsilon         = 10;             % (-) AKAP-mediated scaling factor

%PKAIItot = 0.059;          % (uM) total type 2 PKA % MOUSE
PKAIItot = 0.084; % RABBIT

PKI = PKItot - PKACI_PKI - PKACII_PKI;

dRC_I = - kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I;
dRCcAMP_I = - kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I;
dRCcAMPcAMP_I = - kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI;
dRcAMPcAMP_I = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I;
dPKACI = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI; % ydot(17)
dPKACI_PKI = - kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI;

dRC_II = - kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II;
dRCcAMP_II = - kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II;
dRCcAMPcAMP_II = - kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII;
dRcAMPcAMP_II = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II;
dPKACII = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI; % ydot(18)
dPKACII_PKI = - kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI;
% end PKA module
%% I-1/PP1 module

PP1tot = pin(2); % PP1tot = 0.89; % (uM) total phosphatase 1
I1tot           = 0.3;            % (uM) total inhibitor 1
k_PKA_I1        = 60e-3;          % (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
Km_PKA_I1       = 1.0;            % (uM) Km for I-1 phosphorylation by type 1 PKA
Vmax_PP2A_I1    = 14.0e-3;        % (uM/ms) Vmax for I-1 dephosphorylation by PP2A
Km_PP2A_I1      = 1.0;            % (uM) Km for I-1 dephosphorylation by PP2A

Ki_PP1_I1       = 1.0e-3;         % (uM) Ki for PP1 inhibition by I-1
kf_PP1_I1       = 1;              % (uM) Ki for PP1 inhibition by I-1
kr_PP1_I1       = Ki_PP1_I1;      % (uM) Ki for PP1 inhibition by I-1

I1 = I1tot - I1ptot;
PP1 = PP1tot - I1p_PP1;
I1p = I1ptot - I1p_PP1;
I1_phosph = k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1);
I1_dephosph = Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot);

dI1p_PP1 = kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1;
dI1ptot = I1_phosph - I1_dephosph; % ydot
% end I-1/PP1 module
%% PLB module

PLBtot = pin(3); %p(41) = PLBtot; % PLBtot    [uM]
k_PKA_PLB = 54e-3;     %p(44) = 54;     % k_pka_plb     [1/ms]
Km_PKA_PLB = 21;     %p(45) = 21;     % Km_pka_plb    [uM]
k_PP1_PLB = 8.5e-3;    %p(46) = 8.5;    % k_pp1_plb     [1/ms]
Km_PP1_PLB = 7.0;    %p(47) = 7.0;    % Km_pp1_plb    [uM]

PLB = PLBtot - PLBp;
PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB);
PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp);
dPLBp = PLB_phosph - PLB_dephosph; % ydot
% end PLB module
%% PLM module (from PLB, different total concentration)

PLMtot = pin(4); % p(102) = PLMtot; % PLMtot    [uM]
k_PKA_PLM = 54e-3; % p(103) = 54;     % k_pka_plb     [1/ms]
Km_PKA_PLM = 21; % p(104) = 21;     % Km_pka_plb    [uM]
k_PP1_PLM = 8.5e-3; % p(105) = 8.5;    % k_pp1_plb     [1/ms]
Km_PP1_PLM = 7.0; % p(106) = 7.0;    % Km_pp1_plb    [uM]

PLM = PLMtot - PLMp;
PLM_phosph = k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM);
PLM_dephosph = k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp);
dPLMp = PLM_phosph - PLM_dephosph; % ydot
% end PLM module
%% LCC module

LCCtot = pin(5); %p(53) = LCCtot; % LCCtot        [uM]
PKACII_LCCtot = 0.025;  %p(54) = 0.025;  % PKAIIlcctot   [uM]
PP1_LCC = 0.025;  %p(55) = 0.025;  % PP1lcctot     [uM]
PP2A_LCC = 0.025;  %p(56) = 0.025;  % PP2Alcctot    [uM]
k_PKA_LCC = 54e-3;     %p(57) = 54;     % k_pka_lcc     [1/ms]
Km_PKA_LCC = 21;     %p(58) = 21;%*1.6;     % Km_pka_lcc    [uM]
k_PP1_LCC = 8.52e-3;   %p(59) = 8.52;   % k_pp1_lcc     [1/ms] RABBIT, MOUSE
                %p(59) = 8.5;   % k_pp1_lcc     [1/sec] RAT
Km_PP1_LCC = 3;      %p(60) = 3;      % Km_pp1_lcc    [uM]
k_PP2A_LCC = 10.1e-3;   %p(61) = 10.1;   % k_pp2a_lcc    [1/ms]
Km_PP2A_LCC = 3;      %p(62) = 3;      % Km_pp2a_lcc   [uM]

PKACII_LCC = (PKACII_LCCtot/PKAIItot)*PKACII;
LCCa = LCCtot - LCCap;
LCCa_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa);
LCCa_dephosph = epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap);
dLCCap = LCCa_phosph - LCCa_dephosph; % ydot
LCCb = LCCtot - LCCbp;
LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb);
LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp);
dLCCbp = LCCb_phosph - LCCb_dephosph; % ydot
% end LCC module
%% RyR module

RyRtot = pin(6); %p(63) = RyRtot; % RyRtot        [uM]
PKAIIryrtot = 0.034;  %p(64) = 0.034;  % PKAIIryrtot   [uM]
PP1ryr = 0.034;  %p(65) = 0.034;  % PP1ryr        [uM]
PP2Aryr = 0.034;  %p(66) = 0.034;  % PP2Aryr       [uM]
kcat_pka_ryr = 54e-3;     %p(67) = 54;     % kcat_pka_ryr  [1/ms]
Km_pka_ryr = 21;     %p(68) = 21;     % Km_pka_ryr    [uM]
kcat_pp1_ryr = 8.52e-3;   %p(69) = 8.52;   % kcat_pp1_ryr  [1/ms]
Km_pp1_ryr = 7;      %p(70) = 7;      % Km_pp1_ryr    [uM]
kcat_pp2a_ryr = 10.1e-3;   %p(71) = 10.1;   % kcat_pp2a_ryr [1/ms]
Km_pp2a_ryr = 4.1;    %p(72) = 4.1;    % Km_pp2a_ryr   [uM]

PKACryr = (PKAIIryrtot/PKAIItot)*PKACII;
RyR = RyRtot-RyRp;
RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyR/(Km_pka_ryr+epsilon*RyR);
RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyRp/(Km_pp1_ryr+epsilon*RyRp);
RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyRp/(Km_pp2a_ryr+epsilon*RyRp);
dRyRp = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A; % ydot
% end RyR module
%% TnI module

TnItot = pin(7); %p(73) = TnItot; % TnItot        [uM]
PP2A_TnI = 0.67;   % PP2Atni       [uM]
k_PKA_TnI = 54e-3;     % kcat_pka_tni  [1/ms]
Km_PKA_TnI = 21;     % Km_pka_tni    [uM]
k_PP2A_TnI = 10.1e-3;   % kcat_pp2a_tni [1/ms]
Km_PP2A_TnI = 4.1;    % Km_pp2a_tni   [uM]

TnIn = TnItot - TnIp;
TnI_phosph = k_PKA_TnI*PKACI*TnIn/(Km_PKA_TnI+TnIn);
TnI_dephosph = k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp);
dTnIp = TnI_phosph - TnI_dephosph; % ydot(26)
% end TnI module
%% Myofilament module (from TnI)

Myotot = pin(8);        % MYOtot_bar    [uM]
PP2A_myo = 0.67;             % PP2Amyo       [uM]
kcat_pka_myo = 54e-3;          % kcat_pka_myo  [1/ms]
Km_pka_myo = 21;            % Km_pka_myo    [uM]
kcat_pp2a_myo = 10.1e-3;       % kcat_pp2a_myo [1/ms]
Km_pp2a_myo = 4.1;          % Km_pp2a_myo   [uM]

Myon = Myotot-Myop;  % Non-phos = tot - phos
MyoPHOSPH = kcat_pka_myo*PKACI*Myon/(Km_pka_myo+Myon);
MyoDEPHOSPH = kcat_pp2a_myo*PP2A_myo*Myop/(Km_pp2a_myo+Myop);
dMyop = MyoPHOSPH-MyoDEPHOSPH; % ydot
% end Myo module
%% IKs module

Kstot = pin(9);        % IKstot_bar    [uM]
Yotiao_tot = 0.025;         % Yotiao_tot    [uM]
K_yotiao = 0.1e-3;          % K_yotiao      [uM] ** apply G589D mutation here: 1e2 **
PKAII_ikstot = 0.025;       % PKAII_ikstot  [uM]
PP1_ikstot = 0.025;         % PP1_ikstot    [uM]
k_pka_iks = 1.87e-3; %54;      % k_pka_iks     [1/ms] % adjusted as in Xie et al 2013
Km_pka_iks = 21;            % Km_pka_iks    [uM]
k_pp1_iks = 0.19e-3; %8.52;    % k_pp1_iks     [1/ms] % adjusted as in Xie et al 2013
Km_pp1_iks = 7;             % Km_pp1_iks    [uM]

% Effect of G589D mutation (IKs-Yotiao)
y1 = (-(K_yotiao+Kstot-Yotiao_tot)+sqrt((K_yotiao+Kstot-Yotiao_tot)^2+4*K_yotiao*Yotiao_tot))/2;
x1 = Kstot/(1+y1/K_yotiao);
y2 = (-(K_yotiao+Kstot-Yotiao_tot)-sqrt((K_yotiao+Kstot-Yotiao_tot)^2+4*K_yotiao*Yotiao_tot))/2;
x2 = Kstot/(1+y2/K_yotiao);

free_IKs = x1*(y1 > 0) + x2*(y1 <= 0);
free_Yotiao = y1*(y1 > 0) + y2*(y1 <= 0);
IksYot = free_IKs*free_Yotiao/K_yotiao; % [uM] % IKs-Yotiao

PKACiks = IksYot/Kstot*(PKAII_ikstot/PKAIItot)*PKACII;
PP1iks = IksYot/Kstot*PP1_ikstot;

KSn = Kstot-KSp;  % Non-phos = tot - phos
IKS_PHOSPH = epsilon*k_pka_iks*PKACiks*KSn/(Km_pka_iks+epsilon*KSn);
IKS_DEPHOSPH = epsilon*k_pp1_iks*PP1iks*KSp/(Km_pp1_iks+epsilon*KSp);
dKSp = IKS_PHOSPH-IKS_DEPHOSPH; % ydot
% end Iks module
%% IKr module (from IKs, without mutation)

Krtot = pin(10);        % IKrtot_bar    [uM]
PKAII_ikrtot = 0.025;       % PKAII_ikrtot  [uM]
PP1_ikrtot = 0.025;         % PP1_ikrtot    [uM]
k_pka_ikr = 1.87e-3; %54;      % k_pka_ikr     [1/ms] % adjusted as in Xie et al 2013
Km_pka_ikr = 21;            % Km_pka_ikr    [uM]
k_pp1_ikr = 0.19e-3; %8.52;    % k_pp1_ikr     [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ikr = 7;             % Km_pp1_ikr    [uM]

KRn = Krtot-KRp;  % Non-phos = tot - phos
PKACikr = (PKAII_ikrtot/PKAIItot)*PKACII;
IKR_PHOSPH = epsilon*k_pka_ikr*PKACikr*KRn/(Km_pka_ikr+epsilon*KRn);
IKR_DEPHOSPH = epsilon*k_pp1_ikr*PP1_ikrtot*KRp/(Km_pp1_ikr+epsilon*KRp);
dKRp = IKR_PHOSPH-IKR_DEPHOSPH;
% end IKr module
%% ICFTR module

CFTRtot = pin(11); % ICFTRtot_bar  [uM]
PKAII_CFTRtot = 0.025;      % PKAII_CFTRtot [uM]
PP1_CFTRtot = 0.025;        % PP1_CFTRtot   [uM]
k_pka_CFTR = 54e-3;            % k_pka_CFTR    [1/ms]
Km_pka_CFTR = 8.5;          % Km_pka_CFTR   [uM]
k_pp1_CFTR = 8.52e-3;          % k_pp1_CFTR    [1/ms]
Km_pp1_CFTR = 7;            % Km_pp1_CFTR   [uM]

CFTRn = CFTRtot - CFTRp;  % Non-phos = tot - phos
PKAC_CFTR = (PKAII_CFTRtot/PKAIItot)*PKACII;    % (PKACFTRtot/PKAIItot)*PKAIIact
CFTRphos = epsilon*CFTRn*PKAC_CFTR*k_pka_CFTR/(Km_pka_CFTR+epsilon*CFTRn);
CFTRdephos = PP1_CFTRtot*k_pp1_CFTR*epsilon*CFTRp/(Km_pp1_CFTR+epsilon*CFTRp);
dCFTRp = CFTRphos - CFTRdephos;
% end ICFTR module
%% IClCa module (from CFTR)

ClCatot = pin(12);	% IClCatot_bar  [uM]
PKAII_ClCatot = 0.025;      % PKAII_ClCatot [uM]
PP1_ClCatot = 0.025;        % PP1_ClCatot   [uM]
k_pka_ClCa = 54e-3;            % k_pka_ClCa    [1/ms]
Km_pka_ClCa = 8.5;          % Km_pka_ClCa   [uM]
k_pp1_ClCa = 8.52e-3;          % k_pp1_ClCa    [1/ms]
Km_pp1_ClCa = 7;            % Km_pp1_ClCa   [uM]

ClCan = ClCatot - ClCap;  % Non-phos = tot - phos
PKAC_ClCa = (PKAII_ClCatot/PKAIItot)*PKACII;    % (PKACFTRtot/PKAIItot)*PKAIIact
ClCaphos = epsilon*ClCan*PKAC_ClCa*k_pka_ClCa/(Km_pka_ClCa+epsilon*ClCan);
ClCadephos = PP1_ClCatot*k_pp1_ClCa*epsilon*ClCap/(Km_pp1_ClCa+epsilon*ClCap);
dClCap = ClCaphos - ClCadephos;
% end ICl(Ca) module
%% Ito module (from IKr)

Totot = pin(13);        % Itotot_bar    [uM]
PKAII_itotot = 0.025;       % PKAII_itotot  [uM]
PP1_itotot = 0.025;         % PP1_itotot    [uM]
k_pka_ito = 1.87e-3; %54;      % k_pka_ito     [1/ms] % adjusted as in Xie et al 2013
Km_pka_ito = 21;            % Km_pka_ito    [uM]
k_pp1_ito = 0.19e-3; %8.52;    % k_pp1_ito     [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ito = 7;             % Km_pp1_ito    [uM]

Ton = Totot-Top;  % Non-phos = tot - phos
PKACito = (PKAII_itotot/PKAIItot)*PKACII;
Ito_PHOSPH = epsilon*k_pka_ito*PKACito*Ton/(Km_pka_ito+epsilon*Ton);
Ito_DEPHOSPH = epsilon*k_pp1_ito*PP1_itotot*Top/(Km_pp1_ito+epsilon*Top);
dTop = Ito_PHOSPH-Ito_DEPHOSPH;
% end Ito module
%% IK1 module (from IKr)

K1tot = pin(14);        % IK1tot_bar    [uM]
PKAII_ik1tot = 0.025;       % PKAII_ik1tot  [uM]
PP1_ik1tot = 0.025;         % PP1_ik1tot    [uM]
k_pka_ik1 = 1.87e-3; %54;      % k_pka_ik1     [1/ms] % adjusted as in Xie et al 2013
Km_pka_ik1 = 21;            % Km_pka_ik1    [uM]
k_pp1_ik1 = 0.19e-3; %8.52;    % k_pp1_ik1     [1/ms] % adjusted as in Xie et al 2013
Km_pp1_ik1 = 7;             % Km_pp1_ik1    [uM]

K1n = K1tot-K1p;  % Non-phos = tot - phos
PKACik1 = (PKAII_ik1tot/PKAIItot)*PKACII;
IK1_PHOSPH = epsilon*k_pka_ik1*PKACik1*K1n/(Km_pka_ik1+epsilon*K1n);
IK1_DEPHOSPH = epsilon*k_pp1_ik1*PP1_ik1tot*K1p/(Km_pp1_ik1+epsilon*K1p);
dK1p = IK1_PHOSPH-IK1_DEPHOSPH;
% end IK1 module
%% INa module (from IKr)

Natot = pin(15);        % INatot_bar    [uM]
PKAII_iNatot = 0.025;       % PKAII_iNatot  [uM]
PP1_iNatot = 0.025;         % PP1_iNatot    [uM]
k_pka_iNa = 1.87e-3; %54;      % k_pka_iNa     [1/ms] % adjusted as in Xie et al 2013
Km_pka_iNa = 21;            % Km_pka_iNa    [uM]
k_pp1_iNa = 0.19e-3; %8.52;    % k_pp1_iNa     [1/ms] % adjusted as in Xie et al 2013
Km_pp1_iNa = 7;             % Km_pp1_iNa    [uM]

Nan = Natot-Nap;  % Non-phos = tot - phos
PKACiNa = (PKAII_iNatot/PKAIItot)*PKACII;
INa_PHOSPH = epsilon*k_pka_iNa*PKACiNa*Nan/(Km_pka_iNa+epsilon*Nan);
INa_DEPHOSPH = epsilon*k_pp1_iNa*PP1_iNatot*Nap/(Km_pp1_iNa+epsilon*Nap);
dNap = INa_PHOSPH-INa_DEPHOSPH;
% end INa module
%% ydot

ydot(1)=dLR;
ydot(2)=dLRG;
ydot(3)=dRG;
ydot(4)=db1AR_S464;
ydot(5)=db1AR_S301;
ydot(6)=dGsaGTPtot;
ydot(7)=dGsaGDP;
ydot(8)=dGsby;
ydot(9)=dAC_GsaGTP;
%ydot(10)=dPDEp;
ydot(10)=dPDE3p;
ydot(11)=dPDE4p;
ydot(12)=dcAMPtot;
ydot(13)=dRC_I;
ydot(14)=dRCcAMP_I;
ydot(15)=dRCcAMPcAMP_I;
ydot(16)=dRcAMPcAMP_I;
ydot(17)=dPKACI;
ydot(18)=dPKACI_PKI;
ydot(19)=dRC_II;
ydot(20)=dRCcAMP_II;
ydot(21)=dRCcAMPcAMP_II;
ydot(22)=dRcAMPcAMP_II;
ydot(23)=dPKACII;
ydot(24)=dPKACII_PKI;
ydot(25)=dI1p_PP1; % output CaMKII
ydot(26)=dI1ptot;
ydot(27)=dPLBp; % output
ydot(28)=dPLMp; % output
ydot(29)=dLCCap; % output
ydot(30)=dLCCbp; % output
ydot(31)=dRyRp; % output
ydot(32)=dTnIp; % output
ydot(33)=dMyop; % output
ydot(34)=dKSp; % output
ydot(35)=dKRp; % output
ydot(36)=dCFTRp; % output
ydot(37)=dClCap; % output
ydot(38)=dTop; % output
ydot(39)=dK1p; % output
ydot(40)=dNap; % output
