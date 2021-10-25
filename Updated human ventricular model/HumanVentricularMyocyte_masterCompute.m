%% Main file: HumanVentricularMyocyte_masterCompute
% This file loads the initial conditions, calls the ode solver, and plots
% the results.

close all
clear 
clc
%% Parameters

% CaMKII
camkii_exp = 1; % 0, 1 or 6
% ISO
Ligtot = 0; % 0.1, 0.05, 0.02 or 0 % [uM] SET ISO CONCENTRATION HERE
% Protocol
prot_index = 1; % 0 no stim, 1 pace, 2 v-step
prot_freq = 1; % [Hz] CHANGE DEPENDING ON FREQUENCY
prot_cycleLength = 1e3/prot_freq;     % [ms]
prot_input_par = 10; % Input parameter for stimulation protocols
% Ca clamp
Ca_clamp_index = 0; % 0 Ca-free, 1 Ca-clamp to initial value
% Na clamp
Na_clamp_index = 0; % 0 Na-free, 1 Na-clamp to initial value
% Simulation
flag_ECC = 1;    % if 0, module clamped
flag_cam = 1;    % if 0, module clamped
flag_CaMKII = 1; % if 0, module clamped
flag_BAR = 1;    % if 0, module clamped
% State variables
ny_ECC = 83; % state vars in ECC module
ny_cam = 15; % state vars in each CaM module
ny_CaMKII = 6; % state vars in CaMKII module
ny_BAR = 39; % state vars in BAR module
% Sensitivity Analysis
% 1) GNa 2) GNaL 3) GNaB 4) vNKA 5) Gtof
% 6) Gtos 7) GKs 8) GKr 9) GKur1 10) GKur2
% 11) Gss 12) GKp 13) GK1 14) GCFTR 15) GClCa
% 16) GClB 17) GCaL 18) GCaB 19) vNCX 20) vPMCA
% 21) vSERCA 22) vRel 23) vLeak
par_SA = ones(1,23);

p = [camkii_exp Ligtot prot_index prot_cycleLength prot_input_par...
    Ca_clamp_index Na_clamp_index flag_ECC flag_cam flag_CaMKII...
    flag_BAR ny_ECC ny_cam ny_CaMKII ny_BAR par_SA];
%% Load initial conditions

% NEW - Myofilament unplugged
load yf_hvm_1Hz_NEW
%load yf_hvm_1Hz_120s_ISO_NEW

y0n = yfinal;
%% Define [Ca] for simulations with fixed [Ca]

if Ca_clamp_index == 1
    Ca_conc_clamp = 500*1e-6; % (mM)
    y0n(38) = Ca_conc_clamp;
    y0n(37) = Ca_conc_clamp;
    y0n(36) = Ca_conc_clamp;
end
%% Define [Na] for simulations with fixed [Na], if different from ICs

if Na_clamp_index == 1
    Na_conc_clamp = 10; % (mM)
    y0n(34) = Na_conc_clamp;
end
%% Establish and define globals

global tStep tArray ICa_store Ito_store INa_store IK1_store 
global Jserca_store IKs_store Jleak_store Incx_store
global INaK_store Ikr_store INabk_store
global Lmyo_store Fmyo_store Vmax_store
global IKp_store Iclca_store Iclbk_store Ipca_store Icabk_store Cai_store

tStep = 1; tArray = zeros(1,1e6); ICa_store = zeros(1,1e6); Ito_store = zeros(3,1e6);
INa_store = zeros(1,1e6); IK1_store = zeros(1,1e6); 
Jserca_store = zeros(1,1e6); IKs_store = zeros(1,1e6);
Jleak_store = zeros(1e6,2); Incx_store = zeros(1,1e6);
INaK_store = zeros(1,1e6); Ikr_store = zeros(1,1e6); INabk_store = zeros(1,1e6);
Fmyo_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6); Vmax_store = zeros(1,1e6);
IKp_store = zeros(1,1e6); Iclca_store = zeros(1,1e6); Iclbk_store = zeros(1,1e6);
Ipca_store = zeros(1,1e6); Icabk_store = zeros(1,1e6); Cai_store = zeros(1,1e6);
%% Run single simulation

tic
tspan = [0 3e3]; % [ms]
options = odeset('RelTol',1e-5,'MaxStep',2); 
[t,y] = ode15s(@HumanVentricularMyocyte_masterODEfile,tspan,y0n,options,p);
yfinal = y(end,:)';
toc
%% Saving final conditions

% NEW - Myofilament unplugged
%save yf_hvm_1Hz_NEW yfinal
%save yf_hvm_1Hz_120s_ISO_NEW yfinal
%% Rename outputs

tA = tArray(1:tStep); dVm = Vmax_store(1:tStep);
ICa = ICa_store(1:tStep); Ito = Ito_store(1,1:tStep);
Itof = Ito_store(2,1:tStep); Itos = Ito_store(3,1:tStep);
INa = INa_store(1:tStep); IK1 = IK1_store(1:tStep);
IKs = IKs_store(1:tStep); IKr = Ikr_store(1:tStep); 
IKp = IKp_store(1:tStep);
IClCa = Iclca_store(1:tStep); IClbk = Iclbk_store(1:tStep);
INCX = Incx_store(1:tStep); INaK = INaK_store(1:tStep);
INabk = INabk_store(1:tStep); IPMCA = Ipca_store(1:tStep);
ICabk = Icabk_store(1:tStep); Cai_tA = Cai_store(1:tStep);
Jserca = Jserca_store(1:tStep); Jleak = Jleak_store(1:tStep,:);
Lm = Lmyo_store(1:tStep); Fm = Fmyo_store(1:tStep);
%% Plot
%t = t-600e3; tA = tA-600e3;

figure,set(gcf,'color','w')
subplot(4,2,1); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,3); hold on, plot(t,y(:,38)*1e6); ylabel('[Ca]i (nM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,5); hold on, plot(t,y(:,31)); ylabel('[Ca]SR (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,7); hold on, plot(t,y(:,34)); ylabel('[Na]i (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')

subplot(4,2,2); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,4); hold on, plot(t,y(:,38)*1e6); ylabel('[Ca]i (nM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,6); hold on, plot(tA,ICa); ylabel('ICaL (A/F)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,8); hold on,
plot(tA,Lm); ylabel('Length (um)');
%plot(tA,Fm); ylabel('Force (mN/mm2)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')

figure,set(gcf,'color','w')
subplot(4,2,1); plot(tA,INa); title('INa');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,3); plot(tA,ICa); title('ICa');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,5); plot(tA,INCX); title('INCX');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,7); plot(tA,INaK); title('INaK');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,2); plot(tA,IKr); title('IKr');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,4); plot(tA,IKs); title('IKs');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,6); plot(tA,Ito); title('Ito');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,8); plot(tA,IK1); title('IK1');
set(gca,'box','off','tickdir','out','fontsize',12)

% figure,set(gcf,'color','w')
% hold on, plot(tA,INa); ylabel('INa (A/F)');
% set(gca,'box','off','tickdir','out','fontsize',12)
% xlabel('Time (ms)')

% CaMKII
LCCtotDyad = 31.4*.9;       % [uM] - Total Dyadic [LCC] - (umol/l dyad)
RyRtot = 382.6;             % [uM] - Total RyR (in Dyad)
PLBtot = 38;                % [uM] - Total [PLB] in cytosolic units
LCC_CKp = y(:,ny_ECC+3*ny_cam+2)./LCCtotDyad;   % fractional CaMKII-dependent LCC dyad phosphorylation
RyR_CKp = y(:,ny_ECC+3*ny_cam+4)./RyRtot;       % fractional CaMKII-dependent RyR phosphorylation
PLB_CKp = y(:,ny_ECC+3*ny_cam+5)./PLBtot;       % fractional CaMKII-dependent PLB phosphorylation

figure,set(gcf,'color','w')
subplot(4,1,1); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
set(gca,'box','off','tickdir','out','fontsize',12)
title('CaMKII')
subplot(4,1,2); hold on, plot(t,LCC_CKp*100); ylabel('LTCC (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on, plot(t,RyR_CKp*100); ylabel('RyR/Nav (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on, plot(t,PLB_CKp*100); ylabel('PLB (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')

% bAR
cAMPtot = y(:,ny_ECC+3*ny_cam+ny_CaMKII+12);

LCCtotBA = 0.025;           % [uM] - [umol/L cytosol]
RyRtotBA = 0.135;           % [uM] - [umol/L cytosol]
PLBtotBA = PLBtot;          % [uM] - [umol/L cytosol]
TnItotBA = 70;              % [uM] - [umol/L cytosol]
IKstotBA = 0.025;           % [uM] - [umol/L cytosol]
ICFTRtotBA = 0.025;         % [uM] - [umol/L cytosol]
PLMtotBA = 48;              % [uM] - [umol/L cytosol] as in Yang & Saucerman (mouse) model
MyototBA = 70;              % [uM] - [umol/L cytosol] as TnI
IKrtotBA = 0.025;           % [uM] - [umol/L cytosol] as IKs
IClCatotBA = 0.025;         % [uM] - [umol/L cytosol] as ICFTR
ItototBA = 0.025;         % [uM] - [umol/L cytosol] as IKr
IK1totBA = 0.025;         % [uM] - [umol/L cytosol] as IKr
INatotBA = 0.025;         % [uM] - [umol/L cytosol] as IKr
PLB_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+27)./PLBtotBA;
PLM_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+28)./PLMtotBA;
LCCa_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+29)./LCCtotBA;
LCCb_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+30)./LCCtotBA;
RyR_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+31)./RyRtotBA;
TnI_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+32)./TnItotBA;
Myo_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+33)./MyototBA;
IKs_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+34)./IKstotBA;
IKr_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+35)./IKrtotBA;
IClCa_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+36)./IClCatotBA;
Ito_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+37)./ItototBA;
IK1_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+38)./IK1totBA;
INa_PKAp = y(:,ny_ECC+3*ny_cam+ny_CaMKII+39)./INatotBA;

figure,set(gcf,'color','w')
subplot(4,2,1); hold on, plot(t,cAMPtot); ylabel('cAMP (uM)');
set(gca,'box','off','tickdir','out','fontsize',12)
title('bAR')
subplot(4,2,2); hold on, plot(t,LCCa_PKAp*100,t,LCCb_PKAp*100); ylabel('LTCC (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('a','b')
subplot(4,2,3); hold on, plot(t,RyR_PKAp*100); ylabel('RyR (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,4); hold on, plot(t,PLB_PKAp*100); ylabel('PLB (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,5); hold on, plot(t,PLM_PKAp*100); ylabel('PLM (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,6); hold on, plot(t,Myo_PKAp*100); ylabel('TnI/Myo (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,2,7); hold on, plot(t,IKs_PKAp*100,t,IKr_PKAp*100); ylabel('IKs/IKr (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
legend('IKs','IKr')
xlabel('Time (ms)')
subplot(4,2,8); hold on, plot(t,IClCa_PKAp*100); ylabel('ClCa (%)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')
