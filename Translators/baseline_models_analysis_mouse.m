%% Main file
% This file loads the initial conditions, calls the ode solver, and plots
% the results.

%close all
%clear 
%clc
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
flag_BAR = 0;    % if 0, module clamped
% State variables
ny_ECC = 93; % state vars in ECC module
ny_cam = 15; % state vars in each CaM module
ny_CaMKII = 6; % state vars in CaMKII module
ny_BAR = 36; % state vars in BAR module
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
load yf_mvm_1Hz_NEW
    
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
%% Establish and define globals variables

global tStep tArray ICa_store Ito_store INa_store IK1_store 
global Jserca_store IKs_store Jleak_store ICFTR_store Incx_store
global Iss_store dVm_store Ipca_store INaK_store INabk_store Ikr_store IKp_store
global Ikur1_store Ikur2_store Iclca_store Iclbk_store Icabk_store
global Lmyo_store Fmyo_store Cai_store
 
tStep = 1; tArray = zeros(1,1e6); ICa_store = zeros(1,1e6);
Ito_store = zeros(3,1e6); INa_store = zeros(1,1e6);
IK1_store = zeros(1,1e6); Jserca_store = zeros(1,1e6);
IKs_store = zeros(1,1e6); Jleak_store = zeros(1e6,2);
ICFTR_store = zeros(1,1e6); Incx_store = zeros(1,1e6);
Ikur1_store = zeros(1,1e6); Ikur2_store = zeros(1,1e6);
Iss_store = zeros(1,1e6); dVm_store = zeros(1,1e6);
Ipca_store = zeros(1,1e6); INaK_store = zeros(1,1e6);
INabk_store = zeros(1,1e6); Ikr_store = zeros(1,1e6);
IKp_store = zeros(1,1e6); Iclca_store = zeros(1,1e6);
Iclbk_store = zeros(1,1e6); Icabk_store = zeros(1,1e6);
Cai_store = zeros(1,1e6); Lmyo_store = zeros(1,1e6);
Fmyo_store = zeros(1,1e6);
%% Run simulation

tic
tspan = [0 3e3]; % [ms]
options = odeset('RelTol',1e-5,'MaxStep',2);
[t,y] = ode15s(@MouseVentricularMyocyte_masterODEfile,tspan,y0n,options,p);
y0n_baseline = y(end,:)';
toc
%% Plot

figure,set(gcf,'color','w')
subplot(4,1,1); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
title('Baseline Mouse model (control, 1 Hz)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on, plot(t,y(:,38)*1e6); ylabel('[Ca]i (nM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on, plot(t,y(:,31)); ylabel('[Ca]SR (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on, plot(t,y(:,34)); ylabel('[Na]i (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')

%% Analysis

time = t; % (ms)
Vm = y(:,39); % (mV)
Ca = y(:,38); % (mM)
CaSR = y(:,31); % (mM)
Na = y(:,34); % (mM)
% dVm = currents(:,1); % (mV/ms)
dVm_array = (y(2:end,39)-y(1:end-1,39))./(t(2:end)-t(1:end-1));
dVm = [dVm_array(1); dVm_array];
period = prot_cycleLength; % ms
AP_index = 2; % w/ 1 first AP, otherwise last AP
Ca_clamp = Ca_clamp_index;

% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp
outputs = function_beat_analysis(time,Vm,Ca,CaSR,Na,dVm,period,AP_index,Ca_clamp);

outputs_mouse_baseline = outputs;
%% Ion current block (50%)

parameter_names = {'GNa' 'GNaL' 'GNaB' 'vNKA' 'Gtof'...
    'Gtos' 'GKs' 'GKr' 'GKur1' 'GKur2'...
    'Gss' 'GKp' 'GK1' 'GCFTR' 'GClCa'...
    'GClB' 'GCaL' 'GCaB' 'vNCX' 'vPMCA'...
    'vSERCA' 'vRel' 'vLeak'} ;
% 1) GNa 2) GNaL 3) GNaB 4) vNKA 5) Gtof
% 6) Gtos 7) GKs 8) GKr 9) GKur1 10) GKur2
% 11) Gss 12) GKp 13) GK1 14) GCFTR 15) GClCa
% 16) GClB 17) GCaL 18) GCaB 19) vNCX 20) vPMCA
% 21) vSERCA 22) vRel 23) vLeak

% Currents/Transporters to block
block_index = [2 5 8 13 17 19 21 22];
block_factor = 0.5;

yfinal_matrix = zeros(length(block_index),length(y0n_baseline));
output_matrix = zeros(length(block_index),18);

tic
parfor ii = 1:length(block_index)
    par_SA = ones(1,23);
    par_block_index = block_index(ii);
    par_SA(par_block_index) = block_factor;
    
    p = [camkii_exp Ligtot prot_index prot_cycleLength prot_input_par...
        Ca_clamp_index Na_clamp_index flag_ECC flag_cam flag_CaMKII...
    	flag_BAR ny_ECC ny_cam ny_CaMKII ny_BAR par_SA];
    
    tspan = [0 300e3]; % [ms]
    [t,y] = ode15s(@MouseVentricularMyocyte_masterODEfile,tspan,y0n_baseline,options,p);
    yfinal_matrix(ii,:) = y(end,:)';
    
%     figure,set(gcf,'color','w')
%     subplot(4,1,1); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
%     set(gca,'box','off','tickdir','out','fontsize',12)
%     subplot(4,1,2); hold on, plot(t,y(:,38)*1e6); ylabel('[Ca]i (nM)');
%     set(gca,'box','off','tickdir','out','fontsize',12)
%     subplot(4,1,3); hold on, plot(t,y(:,31)); ylabel('[Ca]SR (mM)');
%     set(gca,'box','off','tickdir','out','fontsize',12)
%     subplot(4,1,4); hold on, plot(t,y(:,34)); ylabel('[Na]i (mM)');
%     set(gca,'box','off','tickdir','out','fontsize',12)
%     xlabel('Time (ms)')
end
toc
%% Get outputs

tic
for ii = 1:length(block_index)
    par_SA = ones(1,23);
    par_block_index = block_index(ii);
    par_SA(par_block_index) = block_factor;
    
    p = [camkii_exp Ligtot prot_index prot_cycleLength prot_input_par...
        Ca_clamp_index Na_clamp_index flag_ECC flag_cam flag_CaMKII...
    	flag_BAR ny_ECC ny_cam ny_CaMKII ny_BAR par_SA];
    
    tspan = [0 3e3]; % [ms]
    [t,y] = ode15s(@MouseVentricularMyocyte_masterODEfile,tspan,yfinal_matrix(ii,:),options,p);
    
    figure,set(gcf,'color','w')
    subplot(4,1,1); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
    set(gca,'box','off','tickdir','out','fontsize',12)
    title(parameter_names{par_block_index})
    subplot(4,1,2); hold on, plot(t,y(:,38)*1e6); ylabel('[Ca]i (nM)');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,1,3); hold on, plot(t,y(:,31)); ylabel('[Ca]SR (mM)');
    set(gca,'box','off','tickdir','out','fontsize',12)
    subplot(4,1,4); hold on, plot(t,y(:,34)); ylabel('[Na]i (mM)');
    set(gca,'box','off','tickdir','out','fontsize',12)
    xlabel('Time (ms)')
    
    time = t; % (ms)
    Vm = y(:,39); % (mV)
    Ca = y(:,38); % (mM)
    CaSR = y(:,31); % (mM)
    Na = y(:,34); % (mM)
    % dVm = currents(:,1); % (mV/ms)
    dVm_array = (y(2:end,39)-y(1:end-1,39))./(t(2:end)-t(1:end-1));
    dVm = [dVm_array(1); dVm_array];
    period = prot_cycleLength; % ms
    AP_index = 2; % w/ 1 first AP, otherwise last AP
    Ca_clamp = Ca_clamp_index;

    % 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
    % 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
    % 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
    % 16) CaSRmax 17) CaSRmin 18) CaSRamp
    outputs = function_beat_analysis(time,Vm,Ca,CaSR,Na,dVm,period,AP_index,Ca_clamp);
    
    output_matrix(ii,:) = outputs;
end
toc

index_mouse_block = block_index;
outputs_mouse_block = output_matrix;
%% ISO

load yf_mvm_1Hz_120s_ISO_NEW
y0n = yfinal;

% ISO
Ligtot = 0.1; % 0.1, 0.05, 0.02 or 0 % [uM] SET ISO CONCENTRATION HERE
% Simulation
flag_ECC = 1;    % if 0, module clamped
flag_cam = 1;    % if 0, module clamped
flag_CaMKII = 1; % if 0, module clamped
flag_BAR = 1;    % if 0, module clamped
% Sensitivity Analysis
par_SA = ones(1,23);
    
p = [camkii_exp Ligtot prot_index prot_cycleLength prot_input_par...
    Ca_clamp_index Na_clamp_index flag_ECC flag_cam flag_CaMKII...
    flag_BAR ny_ECC ny_cam ny_CaMKII ny_BAR par_SA];

tic
tspan = [0 3e3]; % [ms]
[t,y] = ode15s(@MouseVentricularMyocyte_masterODEfile,tspan,y0n,options,p);
toc

% Plot
figure,set(gcf,'color','w')
subplot(4,1,1); hold on, plot(t,y(:,39)); ylabel('Em (mV)');
title('Baseline Mouse model (ISO, 1 Hz)')
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,2); hold on, plot(t,y(:,38)*1e6); ylabel('[Ca]i (nM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,3); hold on, plot(t,y(:,31)); ylabel('[Ca]SR (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
subplot(4,1,4); hold on, plot(t,y(:,34)); ylabel('[Na]i (mM)');
set(gca,'box','off','tickdir','out','fontsize',12)
xlabel('Time (ms)')

% Analysis
time = t; % (ms)
Vm = y(:,39); % (mV)
Ca = y(:,38); % (mM)
CaSR = y(:,31); % (mM)
Na = y(:,34); % (mM)
% dVm = currents(:,1); % (mV/ms)
dVm_array = (y(2:end,39)-y(1:end-1,39))./(t(2:end)-t(1:end-1));
dVm = [dVm_array(1); dVm_array];
period = prot_cycleLength; % ms
AP_index = 2; % w/ 1 first AP, otherwise last AP
Ca_clamp = Ca_clamp_index;

outputs = function_beat_analysis(time,Vm,Ca,CaSR,Na,dVm,period,AP_index,Ca_clamp);

outputs_mouse_iso = outputs;
%% Saving

%save outputs_current_block_mouse index_mouse_block outputs_mouse_baseline outputs_mouse_block outputs_mouse_iso
