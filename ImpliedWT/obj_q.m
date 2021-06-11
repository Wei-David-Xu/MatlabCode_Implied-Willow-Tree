function MSE = obj_q(q, S, K, call_mkt, put_mkt, B0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Objective function when optimize probabilities of each nodes in willow 
%   tree, i.e., Mean-Square-Error of |V_mkt - V_mdl|
%
% Input:
%   q: probabilities of each nodes in willow tree
%   S: willow tree of stock price
%   K : strike price vector
%   call_mkt: vector of call option prices in the market
%   put_mkt : vector of put  option prices in the market
%   B0: discount factor, i.e., exp(-r*T)
%
% Output:
%   MSE: Mean-Square-Error of |V_mkt - V_mdl|
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K);
SE = 0;
for n_k = 1:Num_K
    SE = SE + (call_mkt(n_k) - B0* q'*max(S - K(n_k), 0))^2;
    SE = SE + (put_mkt(n_k) - B0* q'*max(K(n_k)-S, 0))^2;
end
MSE = SE/Num_K/2;
end
