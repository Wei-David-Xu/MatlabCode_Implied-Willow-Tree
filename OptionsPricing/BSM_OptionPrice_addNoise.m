function [call_mkt, put_mkt] = BSM_OptionPrice_addNoise(S0,K,r,sigma,time_nodes, N, noise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Generate European options prices under GBM by Black-Scholes formula 
%   adding noise.
%
% Input:
%   S0: initial stock price
%   K : strikes (or vector of strikes)
%   r : risk-free interest rate
%   time_nodes : all maturities of options
%   N : the number of time steps
%   sigma : volatility
%
% Output:
%   call_mkt: European call options prices adding noise
%   put_mkt : European put options prices adding noise
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K);
call_mkt = [];
put_mkt = [];
for n_t = 1:N
    T_n = time_nodes(n_t);
    % Black-Scholes formula
    [call, put] = BS_GBM_Euro(S0,K,r,T_n,sigma);
    % Add noise for Options price
    call = call + noise*call .* randn(1, Num_K);
    put = put + noise*put .* randn(1, Num_K);
    call_mkt = [call_mkt, call'];
    put_mkt = [put_mkt, put'];
end
end