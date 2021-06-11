function [S, p, q] = ImpWT_given_moments_underQ(S0, B0, K, call_mkt, put_mkt,...
                         m, N, gamma, mean_dd, var_dd, skew_dd, kurt_dd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Construct implied willow tree under Q-measure given the options prices 
%   and moments under Q-measure.
% 
% Input:
%   S0: initial stock price
%   B0: discount factor, i.e., exp(-r*time_nodes) vector
%   K : strike price vector
%   call_mkt: matrix of call option prices in the market under Q-measure 
%             (number of strikes * number of maturities)
%   put_mkt : matrix of put option prices in the market under Q-measure 
%             (number of strikes * number of maturities)
%   m : the number of willow tree nodes
%   N : the number of time steps
%   gamma: parameter for willow tree
%   mean_dd : mean of log return under Q-measure 
%   var_dd  : variance of log return under Q-measure
%   skew_dd : skewness of log return under Q-measure 
%   kurt_dd : kurtosis of log return under Q-measure 
%
% Output:
%   S : willow tree of stock price under Q-measure
%   p : transition probability matrix under Q-measure
%   q : probabilities of each nodes in willow tree under Q-measure
%
% Implemented by
%      Bing Dong at 2021.04
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate asset price of each nodes S(i,n)
G = [mean_dd; var_dd; skew_dd; kurt_dd];
[~,~,S,~] = WTnodes_from_JohnsonCurve(G,N,m,gamma);
S = S0*exp(S);
S = sort(S);

% Calculate the probabilities of each nodes q(i,n)
q0 = 1/m*ones(m,N);
q = 1/m*ones(m,N);
for n_t = 1:N
    x0 = q0(:, n_t);
    Aeq = [ones(1,m); S(:,n_t)'./S0*B0(n_t)];
    Beq = [1; 1];
    LB = zeros(m,1);
    UB = ones(m,1);
    options = optimoptions(@fmincon, 'Display','iter','Algorithm','sqp',...
                  'MaxFunEval',inf,'MaxIter',Inf);
    [qn_opt, fval] = fmincon(@(qq)obj_q(qq, S(:,n_t), K, ...
        call_mkt(:,n_t), put_mkt(:,n_t), B0(n_t)),...
        x0,[],[],Aeq,Beq,LB,UB,[],options);
    q(:, n_t) = qn_opt;
end

% Calculate transition probabilities p(i,j,n)
p0 = 1/m*ones(m*m,1);
p = zeros(m,m,N-1);
for n_t = 1:N-1
    Aeq_Metrix = [ones(1,m); S(:,n_t+1)'*B0(n_t+1)/B0(n_t)];
    Aeq_log = repmat('Aeq_Metrix,',1,m);
    Aeq = eval(sprintf('blkdiag(%s)',Aeq_log(1:end-1)));
    Beq = [ones(1,m); S(:,n_t)'];
    Beq = reshape(Beq, 2*m,1);
    LB = zeros(m*m, 1);
    UB = ones(m*m, 1);
    [p_opt, fval] = fmincon(@(pp)obj_p_givenq(pp, q(:, n_t),q(:, n_t+1)),...
        p0, [], [], Aeq, Beq, LB, UB, [], options); 
    p_opt = (reshape(p_opt,m,m))';
    p(:,:,n_t) = p_opt;
end
end


