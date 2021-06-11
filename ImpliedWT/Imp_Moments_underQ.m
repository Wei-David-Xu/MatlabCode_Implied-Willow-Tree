function [mean_dd, var_dd, skew_dd, kurt_dd] = Imp_Moments_underQ(call_mkt, put_mkt, S0, B0, K, dK,N) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate the implied moments matrix using out-of-money option price 
%   under Q-measure.
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
%
% Output:
%   mean_dd : mean of log return under Q-measure 
%   var_dd  : variance of log return under Q-measure
%   skew_dd : skewness of log return under Q-measure 
%   kurt_dd : kurtosis of log return under Q-measure 
%
% Implemented by
%      Bing Dong at  April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n_t = 1:N
    % Implied Moments at each time step
    [mean_dd(n_t), var_dd(n_t), skew_dd(n_t), kurt_dd(n_t)] = ...
        Imp_Moments_underQ_oneMaturity(call_mkt(:,n_t)', put_mkt(:,n_t)', S0, B0(n_t), K, dK);
end
end

function [mean_dd, var_dd, skew_dd, kurt_dd] = ...
    Imp_Moments_underQ_oneMaturity(call, put, S0, B0_T, K, dK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate the implied moments matrix using out-of-money option price 
%   for one maturity under Q-measure .
%
% Input:
%   call: vector of call option prices under Q-measure (1 * number of strikes)
%   put : vector of put  option prices under Q-measure (1 * number of strikes)
%   S0: initial stock price
%   B0_T: discount factor , i.e., exp(-r*T)
%   K : strikes (or vector of strikes)
%   dK: delta K
%   N : the number of time steps
%
% Output:
%   mean_dd : mean of log return under Q-measure 
%   var_dd  : variance of log return under Q-measure
%   skew_dd : skewness of log return under Q-measure 
%   kurt_dd : kurtosis of log return under Q-measure 
%
% Implemented by
%      Bing Dong at  April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K0 = S0;
Num_K = length(K);
out_of_money = zeros(1, Num_K);
for i = 1:Num_K
    if K(i) < K0
        out_of_money(i) = put(i);
    else
        out_of_money(i) = call(i);
    end
end

k1 = (log(K0/S0) + S0/(B0_T*K0) - 1) - sum( dK./K.^2./B0_T .* out_of_money);

k2 = (log(K0/S0))^2 + 2*log(K0/S0)/K0*(S0/B0_T-K0) ...
    + sum( dK*2./K.^2 .* (1-log(K./S0))./B0_T .* out_of_money);

k3 = (log(K0/S0))^3 + 3*(log(K0/S0))^2/K0*(S0/B0_T-K0) ...
    + sum( dK*3./K.^2 .*( log(K./K0).*(2-log(K./K0)) ...
    + 2*log(K0/S0).*(1-log(K./K0)) ...
    - (log(K0/S0))^2 )./B0_T .* out_of_money );

k4 = (log(K0/S0))^4 + 4*(log(K0/S0))^3/K0*(S0/B0_T-K0) ...
    + sum( dK*4./K.^2.*( (log(K./K0)).^2.*(3-log(K./K0)) ...
    + 3*log(K0/S0).*log(K./K0).*(2-log(K./K0))...
    + 3*(log(K0/S0)^2.*(1-log(K./K0)))...
    - (log(K0/S0))^3 ) ./B0_T.* out_of_money );

mean_dd = k1;
var_dd = k2 - k1.^2;
skew_dd = (k3 - 3*k1*var_dd - k1^3)/var_dd^(3/2);
kurt_dd = (k4 - 4*k3*k1 + 6*k2*k1^2 - 3*k1^4)/var_dd^2;

end