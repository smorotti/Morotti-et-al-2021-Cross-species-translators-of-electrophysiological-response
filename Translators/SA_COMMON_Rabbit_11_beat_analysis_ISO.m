% This file analyzes the properties of each model in the population.
% It creates the matrix 'all_outputs' and the arrays 'output_names' and
% 'output_units'.

close all
clear
clc

%% Loading initial conditions
% load matrix all_ICs (columns: N state variables, rows: N trials)
load ICs_matrix_1500_0p1_60s_ISO_rabbit_COMMON_3p5Hz
    
% columns: N state variables
% rows: N trials

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load parameter_matrix_1500_COMMON

[N_trials N_par] = size(all_parameters);

%% Simulation parameters
% CaMKII
camkii_exp = 1; % 0, 1 or 6
% ISO
Ligtot = 0.1; % 0.1, 0.05, 0.02 or 0 % [uM] SET ISO CONCENTRATION HERE
% Protocol
prot_index = 1; % 0 no stim, 1 pace, 2 v-step
prot_freq = 3.5; % [Hz] CHANGE DEPENDING ON FREQUENCY
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
ny_BAR = 40; % state vars in BAR module
% Sensitivity Analysis
% 1) GNa 2) GNaL 3) GNaB 4) vNKA 5) Gtof
% 6) Gtos 7) GKs 8) GKr 9) GKur1 10) GKur2
% 11) Gss 12) GKp 13) GK1 14) GCFTR 15) GClCa
% 16) GClB 17) GCaL 18) GCaB 19) vNCX 20) vPMCA
% 21) vSERCA 22) vRel 23) vLeak
par_SA = ones(1,23);

tspan = [0 2e3]; % [ms]
%tspan = [0 3*30*2e3]; % [ms]
options = odeset('RelTol',1e-5,'MaxStep',2);

%% Outputs
% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp
        
output_names = {'UV', 'APpeak', '|MDP|', 'APamp', 'APD90',...
    'APD70', 'APD50', 'APD30', 'CaTmax', 'CaTmin',...
    'CaTamp', 'CaTttp', 'CaTt50', 'CaTtau', 'Na',...
    'CaSRmax', 'CaSRmin', 'CaSRamp'};

output_units = {'mV/ms', 'mV', 'mV', 'mV', 'ms',...
    'ms', 'ms', 'ms', 'mM', 'mM',...
    'mM', 'ms', 'ms', 'ms', 'mM',...
    'mM', 'mM', 'mM'};

N_outputs = length(output_names); % number of outputs of beat analysis

%% Run cycle
all_outputs = zeros(N_trials,N_outputs);

tic
parfor ii = 1:N_trials
%for ii = 1:N_trials % plot figure
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    y0n = all_ICs(ii,:);
    
    par_SA = all_parameters(ii,:); % 23 parameters
    p = [camkii_exp Ligtot prot_index prot_cycleLength prot_input_par...
        Ca_clamp_index Na_clamp_index flag_ECC flag_cam flag_CaMKII...
        flag_BAR ny_ECC ny_cam ny_CaMKII ny_BAR par_SA];
    
    [t,y] = ode15s(@RabbitVentricularMyocyte_masterODEfile,tspan,y0n,options,p);
        
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
    all_outputs(ii,:) = outputs;
end

all_outputs
% columns: N outputs
% rows: N trials
toc

%% Saving
%save outputs_matrix_1500_0p1_60s_ISO_rabbit_COMMON_3p5Hz all_outputs output_names output_units
    