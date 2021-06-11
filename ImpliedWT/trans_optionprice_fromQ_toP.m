function [call_P, put_P] = trans_optionprice_fromQ_toP(gam,S, q, p, N, S0, B0, K, dK, Num_K, call_mkt, put_mkt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate options prices under P-measure according to the options prices
%   under Q-measure given CRRA-utility.
%     -------------------------------------------------------
%     CRRA Utility
%     U(W) = W^(1-gam)/(1-gam)  if gam > 0 & gam != 1
%            \ln W              if gam = 1
%     -------------------------------------------------------
%
% Input:
%   gam : parameter for CRRA-utility (the magnitude of relative risk aversion)
%   S: willow tree of stock price under Q-measure
%   q: probabilities of each nodes in willow tree under Q-measure
%   p: transition probabilities under Q-measure
%   N : the number of time steps
%   S0: initial stock price
%   B0: discount factor, i.e., exp(-r*time_nodes) vector
%   K : strike price vector
%   dK: delta K
%   call_mkt: call options prices in the market under Q-measure
%   put_mkt : put  options prices in the market under Q-measure
%
% Output:
%   call_P: matrix of call options prices under P-measure
%   put_P : matrix of put  options prices under P-measure
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

payoff = 1; % Digital Option
delta_K1 = 0.1;
K_digital = S0*(delta_K1:delta_K1:1e1);
call_P = [];
put_P = [];

% At t(1)
call_Digital_ddWT = B0(1)*q(:,1)'* payoff*(S(:,1)>K_digital);
% c = \int_0^\inf gam x^(gam-1) e^(r \tau) D(t,\tau,x) dx
%   = \sum_{K_i} gam K_i^(gam-1) e^(r \tau) D(t,\tau,K_i)
c = gam * B0(1)^(-1) * K_digital * call_Digital_ddWT' * delta_K1*S0;
call_Digital_ddWT_K = B0(1)*q(:,1)'* payoff*(S(:,1)>K);
for i = 1:Num_K
    K_phy = K(i);
    call_tn_phy(i) = call_mkt(i,1)/c./CRRA_diff1(K_phy, gam) + ...
        CRRA_diff1_2(K,gam)/c .* ( call_mkt(:,1)' + ...
        (K-K_phy).*call_Digital_ddWT_K ) .* dK * (K>K_phy)';
    put_tn_phy(i) = put_mkt(i,1)/c./CRRA_diff1(K_phy, gam) - ...
        CRRA_diff1_2(K,gam)/c .* ( put_mkt(:,1)' + ...
        (K_phy-K).*(B0(1)-call_Digital_ddWT_K) ) .* dK * (K<K_phy)';
end
call_P = [call_P, call_tn_phy'];
put_P = [put_P, put_tn_phy'];

% At t(n_t)
for n_t = 2:N
    [call_Digital_ddWT, ~] = Price_DigitalOption_WT(...
        S(:,1:n_t), p(:,:,1:n_t-1), q(:,1), B0(1:n_t), K_digital, payoff);
    c = gam * B0(n_t)^(-1) * K_digital * call_Digital_ddWT' * delta_K1*S0;
    
    call_Digital_ddWT_K = Price_DigitalOption_WT(...
        S(:,1:n_t), p(:,:,1:n_t-1), q(:,1), B0(1:n_t), K, payoff);
    for i = 1:Num_K
        K_phy = K(i);
        call_tn_phy(i) = call_mkt(i,n_t)/c./CRRA_diff1(K_phy, gam) + ...
            CRRA_diff1_2(K,gam)/c .* ( call_mkt(:,n_t)' + ...
            (K-K_phy).*call_Digital_ddWT_K ) .* dK * (K>K_phy)';
        put_tn_phy(i) = put_mkt(i,n_t)/c./CRRA_diff1(K_phy, gam) - ...
            CRRA_diff1_2(K,gam)/c .* ( put_mkt(:,n_t)' + ...
            (K_phy-K).*(B0(n_t)-call_Digital_ddWT_K) ) .* dK * (K<K_phy)';
    end
    call_P = [call_P, call_tn_phy'];
    put_P = [put_P, put_tn_phy'];
    
end
end

function y1 = CRRA_diff1(W, gam)
% The first order derivative of CRRA-Utility
y1 = W.^(-gam);
end

function y2 = CRRA_diff2(W, gam)
% The second order derivative of CRRA-Utility
y2 = -gam * W.^(-gam-1);
end

function y = CRRA_diff1_2(W, gam)
% output = - CRRA_diff1/CRRA_diff2^2
y = gam*W.^(gam-1);
end

