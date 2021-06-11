% Copyright (c) 2021 Bing Dong, Wei Xu
clear,clc,close all;
% ------------------------------------------------------------------------
% Construct implied willow tree under Q-measure and P-measure
% 
% In this demo, we generate options prices under GBM by Black-Scholes
% formular. You can replace the options prices in the market by replacing 
% [call_mkt, put_mkt] and strike prices K.
% ------------------------------------------------------------------------

% add necessary paths
addpath('../ImpliedWT')
addpath('../OptionsPricing')
addpath('../Utils')

%% Input call and put options prices
% prarmeters for geometric Brownian motion (GBM)
r = 0.05;
sigma = 0.3;
% parameters for options
S0 = 1;
T = 3/12;
N = 3;
dt = T/N;
time_nodes = dt:dt:T; % all maturities
B0 = exp(-r*time_nodes);
delta_K = 0.1;
K = S0*(0.5:delta_K:1.5); % all strike prices
Num_K = length(K);
dK = S0*delta_K*ones(1,Num_K);
% generate options prices under GBM by Black-Scholes formular
noise = 0; % w/o noise
[call_mkt, put_mkt] = BSM_OptionPrice_addNoise(S0,K,r,sigma,time_nodes, N, noise);

%% Construct implied willow tree under risk-neutral (Q) measure
m = 20; % number of nodes for willow tree at each time step
[S, p, q] = Construct_ImpWT_underQ(call_mkt, put_mkt, S0, B0, K, dK, N, m);

%% Construct implied willow tree under physical (P) measure
gam = 2; % CRRA-Utility parameter
[S_P, p_P, q_P] = Construct_ImpWT_underP(gam, S, q, p, N, m, S0, B0,...
    K, dK, call_mkt, put_mkt);
