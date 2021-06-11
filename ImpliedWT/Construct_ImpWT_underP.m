function [S_P, p_P, q_P] = Construct_ImpWT_underP(gam, S, q, p, N, m, S0, ...
    B0, K, dK, call_mkt, put_mkt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Construct implied willow tree given the options prices under P-measure.
%   (set q0 according to Q-measure when optimize probabilities under P-measure)
%
% Input:
%   gam : parameter for CRRA utility
%   S: willow tree of stock price under Q-measure
%   q: probabilities of each nodes in willow tree under Q-measure
%   p: transition probabilities under Q-measure
%   N : the number of time steps
%   m : the number of willow tree nodes
%   S0: initial stock price
%   B0: discount factor, i.e., exp(-r*time_nodes) vector
%   K : strike price vector
%   dK: delta K
%   call_mkt: call option prices in the market under Q-measure
%   put_mkt : put  option prices in the market under Q-measure
%
% Output:
%   S_P : willow tree of stock price under P-measure
%   p_P : transition probability matrix under P-measure
%   q_P : probabilities of each nodes in willow tree under P-measure
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K);
% Calculate options prices under P-measure
[call_P, put_P] = trans_optionprice_fromQ_toP(gam, S, q, p, N, S0, B0, ...
    K, dK, Num_K, call_mkt, put_mkt);

% Calculate Moments under P measure
[mean_dd_P, var_dd_P, skew_dd_P, kurt_dd_P] = ...
    Imp_Moments_underP(call_P, put_P, S0, B0, K, dK, N);

% Construct Willow Tree under P measure
gamma = 0.6; % parameters for willow tree method
[S_P, p_P, q_P] = ImpWT_given_moments_underP(S0, B0, K, call_P, put_P,...
    m, N, gamma, mean_dd_P, var_dd_P, skew_dd_P, kurt_dd_P, q, gam, S);
end