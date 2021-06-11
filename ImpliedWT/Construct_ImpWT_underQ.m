function [S, p, q] = Construct_ImpWT_underQ(call_mkt, put_mkt, S0, B0, K, dK, N, m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Construct implied willow tree given the options prices under Q-measure.
%
% Input:
%   call_mkt: matrix of call option prices in the market under Q-measure 
%             (number of strikes * number of maturities)
%   put_mkt : matrix of put option prices in the market under Q-measure 
%             (number of strikes * number of maturities)
%   S0: initial stock price
%   B0: vector of discount factor , i.e., exp(-r*time_nodes)
%   K : strikes (or vector of strikes)
%   dK: delta K
%   N : the number of time steps
%   m : the number of willow tree nodes
%
% Output:
%   S : willow tree of stock price under Q-measure
%   p : transition probability matrix under Q-measure
%   q : probabilities of each nodes in willow tree under Q-measure
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Moments under Q-measure
[mean_dd, var_dd, skew_dd, kurt_dd] = Imp_Moments_underQ(call_mkt, put_mkt, S0, B0, K, dK, N);

% Construct Willow Tree under Q-measure
gamma = 0.6; % parameter for willow tree method
[S, p, q] = ImpWT_given_moments_underQ(S0, B0, K, call_mkt, put_mkt,...
    m, N, gamma, mean_dd, var_dd, skew_dd, kurt_dd);

end