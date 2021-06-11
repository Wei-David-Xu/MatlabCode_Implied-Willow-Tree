function [call, put] = Price_AmerOption_WT(S, S0, p, q, B0, K_all)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate American call and put options prices by willow tree method.
%
% Input:
%   S: willow tree of stock price -- m*N
%   S0: initial stock price
%   p: transition probabilities -- m*m*N
%   q: transition probability from t0 to t1 -- m*1
%   B0: m*1 : discount factor, i.e., exp(-r*time_nodes) vector
%   K_all: strike prices  -- n_K*1 
%
% Output:
%   call: American call options prices
%   put : American put  options prices
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_K = length(K_all);
for j = 1:Num_K
    K = K_all(j);
    [~, N] = size(S);
    V_call(:,N) = max(S(:,N)-K, 0);
    V_put(:,N) = max(K-S(:,N), 0);
    for n = N-1:-1:1
        V_call(:,n) = max(max(S(:,n)-K, 0), B0(n+1)/B0(n)*p(:,:,n)*V_call(:,n+1));
        V_put(:,n) = max(max(K-S(:,n), 0), B0(n+1)/B0(n)*p(:,:,n)*V_put(:,n+1));
    end
    call(j) = max(max(S0-K, 0), B0(1)*q'*V_call(:,1));
    put(j) = max(max(K-S0, 0), B0(1)*q'*V_put(:,1));
end
end

