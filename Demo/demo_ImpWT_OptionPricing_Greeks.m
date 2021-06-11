% Copyright (c) 2021 Bing Dong, Wei Xu
clear,clc,close all;
% ------------------------------------------------------------------------
% Calcalate the options prices and greeks calculated by implied willow tree
% including European, American and Asian options.
% ------------------------------------------------------------------------

% add necessary paths
addpath('../ImpliedWT')
addpath('../OptionsPricing')
addpath('../Utils')

%% Input call and put options prices
% prarmeters for GBM
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

%% Options Pricing and Greeks
% Pricing European call and put options
[call_Euro_dd, put_Euro_dd,~,~] = Price_EuroOption_WT(S, p, q(:,1), B0, K);
% Pricing American call and put options
[call_Amer_dd, put_Amer_dd] = Price_AmerOption_WT(S, S0, p, q(:,1), B0, K);
% Pricing Asian call options
[call_Asiancall_dd] = Price_AsianCall_WT(T, r, S, p, q(:,1), S0, K);
% Greeks for European call and put options
[delta_call_dd, gamma_call_dd, theta_call_dd, delta_put_dd, gamma_put_dd, ...
    theta_put_dd] = Greeks_EuroOption(S, p, q(:,1), B0, K, T, S0);



