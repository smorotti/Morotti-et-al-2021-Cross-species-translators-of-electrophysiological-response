% This file generates the random perturbations to model parameters.
% It creates the matrix 'all_parameters' with the perturbations values,
% and the array 'parameter_names' with their definitions (strings).

close all
clear
clc

%% Parameters
% 1) GNa 2) GNaL 3) GNaB 4) vNKA 5) Gtof
% 6) Gtos 7) GKs 8) GKr 9) GKur1 10) GKur2
% 11) Gss 12) GKp 13) GK1 14) GCFTR 15) GClCa
% 16) GClB 17) GCaL 18) GCaB 19) vNCX 20) vPMCA
% 21) vSERCA 22) vRel 23) vLeak
parameter_names = {'GNa' 'GNaL' 'GNaB' 'vNKA' 'Gtof'...
    'Gtos' 'GKs' 'GKr' 'GKur1' 'GKur2'...
    'Gss' 'GKp' 'GK1' 'GCFTR' 'GClCa'...
    'GClB' 'GCaL' 'GCaB' 'vNCX' 'vPMCA'...
    'vSERCA' 'vRel' 'vLeak'} ;

n_parameters = length(parameter_names);
baseline_parameters = ones(1,n_parameters);

%% Random variation
variations = 1500; % number of trials

sigmaG = 0.1*ones(1,n_parameters); % standard deviation for parameters
%sigmaG = 0.2*ones(1,n_parameters); % standard deviation for parameters
%sigmaG = 0.3*ones(1,n_parameters); % standard deviation for parameters

all_parameters = zeros(variations,n_parameters);
for ii = 1:n_parameters
    scaling = exp(sigmaG(ii)*randn(1,variations)) ;
    newparams = baseline_parameters(ii)*scaling ;
    all_parameters(:,ii) = newparams ;
end

%all_parameters % size(all_parameters)
% columns: N parameters
% rows: N trials

%% Non-common parameters
% 6) Gtos 7) GKs 9) GKur1 10) GKur2 % 11) Gss 14) GCFTR
all_parameters(:,6) = ones(1,variations) ;
all_parameters(:,7) = ones(1,variations) ;
all_parameters(:,9) = ones(1,variations) ;
all_parameters(:,10) = ones(1,variations) ;
all_parameters(:,11) = ones(1,variations) ;
all_parameters(:,14) = ones(1,variations) ;

all_parameters % size(all_parameters)
% columns: N parameters
% rows: N trials

%% Saving
%save parameter_matrix_1500_COMMON all_parameters parameter_names
