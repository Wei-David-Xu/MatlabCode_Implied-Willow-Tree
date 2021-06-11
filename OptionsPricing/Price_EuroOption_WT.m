function [call, put, V_call, V_put] = Price_EuroOption_WT(S, p, q, B0, K_all)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate European call and put options prices by willow tree method.
%
% Input:
%   S: willow tree of stock price -- m*N
%   p: transition probabilities -- m*m*N
%   q: transition probability from t0 to t1 -- m*1
%   B0: m*1 : discount factor, i.e., exp(-r*time_nodes) vector
%   K_all: strike prices  -- n_K*1 
%
% Output:
%   call: European call options prices
%   put : European put  options prices
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K_all);
[m, N] = size(S);
V_call = zeros(m,N,Num_K);
V_put = zeros(m,N,Num_K);
for j = 1:Num_K
    K = K_all(j);
    V_call(:,N,j) = max(S(:,N)-K,0);
    V_put(:,N,j) = max(K-S(:,N),0);
    for n = N-1:-1:1
        V_call(:,n,j) = B0(n+1)/B0(n)*p(:,:,n)*V_call(:,n+1,j);
        V_put(:,n,j) = B0(n+1)/B0(n)*p(:,:,n)*V_put(:,n+1,j);
    end
    call(j) = B0(1)*q'*V_call(:,1,j);
    put(j) = B0(1)*q'*V_put(:,1,j);
end
end
