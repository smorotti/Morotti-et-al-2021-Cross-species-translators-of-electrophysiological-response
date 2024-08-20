% This file obtains the initial conditions upon perturbations to model parameters.
% It creates the matrix 'all_ICs' with the state variables in each model in
% the population.

close all
clear
clc

%% Loading initial conditions
load yf_rvm_1Hz_NEW % control
N_state_vars = length(yfinal);

% load matrix all_ICs (columns: N state variables, rows: N trials)
load ICs_matrix_1500_0p1_300s_rabbit_COMMON
all_ICs_1Hz = all_ICs;

% columns: N state variables
% rows: N trials

%% Parameters
% load matrix all_parameters (columns: N parameters, rows: N trials)
load parameter_matrix_1500_COMMON % sigma 0.1

[N_trials N_par] = size(all_parameters);

%% Simulation parameters
% CaMKII
camkii_exp = 1; % 0, 1 or 6
% ISO
Ligtot = 0%; % 0.1, 0.05, 0.02 or 0 % [uM] SET ISO CONCENTRATION HERE
% Protocol
prot_index = 1; % 0 no stim, 1 pace, 2 v-step
prot_freq = 0.5; % [Hz] CHANGE DEPENDING ON FREQUENCY
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
flag_BAR = 1*0%;    % if 0, module clamped
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

tspan = [0 300e3]; % [ms]
options = odeset('RelTol',1e-5,'MaxStep',2);

%% Run cycle
all_ICs = zeros(N_trials,N_state_vars);

tic
parfor ii = 1:N_trials
%for ii = 1:1%N_trials
    X = sprintf('Run %d on %d',ii,N_trials); disp(X)
    y0n = all_ICs_1Hz(ii,:);
    
    par_SA = all_parameters(ii,:); % 23 parameters
    p = [camkii_exp Ligtot prot_index prot_cycleLength prot_input_par...
        Ca_clamp_index Na_clamp_index flag_ECC flag_cam flag_CaMKII...
        flag_BAR ny_ECC ny_cam ny_CaMKII ny_BAR par_SA];
        
    [t,y] = ode15s(@RabbitVentricularMyocyte_masterODEfile,tspan,y0n,options,p);
    all_ICs(ii,:) = y(end,:)';
    
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

all_ICs
% columns: N state variables
% rows: N trials
toc

%% Saving
%save ICs_matrix_1500_0p1_300s_rabbit_COMMON_0p5Hz all_ICs
