function [call, put] = BS_GBM_Euro(S0,K,r,T,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate European call and put options prices by Black-Scholes Formula.
%
% Input:
%   S0: initial stock price
%   K : strikes (or vector of strikes)
%   r : risk-free interest rate
%   T : maturity
%   sigma : volatility
%
% Output:
%   call: European call options prices
%   put : European put options prices
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d1 = (log(S0./K)+(r+1/2*sigma^2)*T)/sigma/sqrt(T);
d2 = d1-sigma*sqrt(T);
call = S0*normcdf(d1) - K.*exp(-r*T).*normcdf(d2);
put = K.*exp(-r*T).*normcdf(-d2) - S0*normcdf(-d1);
end



