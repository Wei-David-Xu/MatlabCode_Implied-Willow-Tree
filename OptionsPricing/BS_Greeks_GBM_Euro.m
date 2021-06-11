function [delta_call, delta_put, theta_call, theta_put, gamma_call, gamma_put] ...
    = BS_Greeks_GBM_Euro(S0,K,r,T,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate Greeks for European call and put options by Black-Scholes Formula.
%
% Input:
%   S0: initial stock price
%   K : strikes (or vector of strikes)
%   r : risk-free interest rate
%   T : maturity
%   sigma : volatility
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

d1 = (log(S0./K)+(r+1/2*sigma^2)*T)/sigma/sqrt(T);
d2 = d1-sigma*sqrt(T);


delta_call = normcdf(d1);
delta_put  = normcdf(d1)-1;

theta_call = -S0*normpdf(d1)*sigma/2/sqrt(T) - r*exp(-r*T)*K.*normcdf(d2);
theta_put  = -S0*normpdf(d1)*sigma/2/sqrt(T) + r*exp(-r*T)*K.*normcdf(-d2);

gamma_call = normpdf(d1)/S0/sigma/sqrt(T);
gamma_put  = normpdf(d1)/S0/sigma/sqrt(T);

delta_call = delta_call';
delta_put  = delta_put'; 
theta_call = theta_call';
theta_put  = theta_put'; 
gamma_call = gamma_call';
gamma_put  = gamma_put'; 

end
