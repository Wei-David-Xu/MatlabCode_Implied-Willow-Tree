function [delta_call, gamma_call, theta_call, delta_put, gamma_put, theta_put] ...
    = Greeks_EuroOption(S, p, q, B0, K, T, S0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate Greeks for European call and put options under willow tree.
%
% Note: If you want to calcaulate the greeks for other types of options,
% you can replace the option price [call_Euro_dd, put_Euro_dd] and value
% matrix [V_call, V_put].
%
%
% Input:
%   S : willow tree of stock price under P-measure
%   p : transition probability matrix under P-measure
%   q : probabilities of each nodes in willow tree under P-measure
%   B0: discount factor, i.e., exp(-r*T)
%   K : strikes (or vector of strikes)
%   T : maturity
%   S0: initial stock price
%
% Output:
%   delta_call: delta of European call option
%   delta_put : delta of European put option
%   theta_call: theta of European call option
%   theta_put : theta of European put option
%   gamma_call: gamma of European call option
%   gamma_put : gamma of European put option
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K);
[call_Euro, put_Euro, V_call, V_put] = Price_EuroOption_WT(S, p, q, B0, K);

% calculate greeks for all strike prices
for j = 1:Num_K
    V_call_Euro = V_call(:,:,j);
    V_put_Euro  = V_put(:,:,j);
    [delta_call(j), gamma_call(j), theta_call(j)] ...
        = Greeks(S,V_call_Euro,S0,call_Euro(j),T);
    [delta_put(j), gamma_put(j), theta_put(j)] ...
        = Greeks(S,V_put_Euro,S0,put_Euro(j),T);
end
end