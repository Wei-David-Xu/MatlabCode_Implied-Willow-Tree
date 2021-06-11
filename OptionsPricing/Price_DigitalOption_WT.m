function [call, put] = Price_DigitalOption_WT(S, p, q, B0, K_all, payoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate Digital call and put options prices by willow tree method.
%
% Input:
%   S: willow tree of stock price -- m*N
%   p: transition probabilities -- m*m*N
%   q: transition probability from t0 to t1 -- m*1
%   B0: m*1 : discount factor, i.e., exp(-r*time_nodes) vector
%   K_all: strike prices  -- n_K*1 
%   payoff:payoff of digital option
%
% Output:
%   call: Digital call options prices
%   put : Digital put  options prices
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K_all);
for j = 1:Num_K
    K = K_all(j);
    [~, N] = size(S);
    V_call(:,N) = payoff*(S(:,N)>K);
    V_put(:,N) = payoff*(S(:,N)<K);
    for n = N-1:-1:1
        V_call(:,n) = B0(n+1)/B0(n)*p(:,:,n)*V_call(:,n+1);
        V_put(:,n) = B0(n+1)/B0(n)*p(:,:,n)*V_put(:,n+1);
    end
    call(j) = B0(1)*q'*V_call(:,1);
    put(j) = B0(1)*q'*V_put(:,1);
end
end

