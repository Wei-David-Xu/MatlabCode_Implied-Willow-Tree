function [S, p, q] = ImpWT_given_moments_underP(S0, B0, K, call_P, put_P,...
    m, N, gamma, mean_dd, var_dd, skew_dd, kurt_dd, q_Q, gam, S_Q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Construct implied willow tree under Q-measure given the options prices 
%   and moments under P-measure.
%
% Input:
%   S0: initial stock price
%   B0: discount factor, i.e., exp(-r*time_nodes) vector
%   K : strike price vector
%   call_P: matrix of call option prices under P-measure 
%           (number of strikes * number of maturities)
%   put_P : matrix of put  option prices under P-measure 
%           (number of strikes * number of maturities)
%   m : the number of willow tree nodes
%   N : the number of time steps
%   gamma: parameter for willow tree
%   mean_dd : mean of log return under P-measure 
%   var_dd  : variance of log return under P-measure
%   skew_dd : skewness of log return under P-measure 
%   kurt_dd : kurtosis of log return under under P-measure 
%   q_Q: probabilities of each nodes in willow tree under Q-measure
%   gam : parameter for CRRA-utility (the magnitude of relative risk aversion)
%   S_Q: willow tree of stock price under Q-measure
%
% Output:
%   S : willow tree of stock price under P-measure
%   p : transition probability matrix under P-measure
%   q : probabilities of each nodes in willow tree under P-measure
%
% Implemented by
%      Bing Dong at 2021.04
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate asset price of each nodes S(i,n)
G = [mean_dd; var_dd; skew_dd; kurt_dd];
[~,~,S,~] = WTnodes_from_JohnsonCurve(G,N,m,gamma);% gamma = 0.6;
S = S0*exp(S);
S = sort(S);

% Set q0 according to Q-measure when optimize probabilities under P-measure
for n=1:N
    q_Q(:,n) = interp1(S_Q(:,n),q_Q(:,n),S(:,n),'spline');
    q_Q(:,n) = max(q_Q(:,n),0);
    q_Q(:,n) = q_Q(:,n)/sum(q_Q(:,n));
end
q0 = q_Q./(S.^(-gam))./repmat(sum(q_Q.*(S.^gam)),m,1);

% Calculate the probabilities of each nodes q(i,n)
q = 1/m*ones(m,N);
for n_t = 1:N
    x0 = q0(:, n_t);
    Aeq = [ones(1,m)];
    Beq = [1];
    LB = zeros(m,1);
    UB = ones(m,1);
    options = optimoptions(...  
          @fmincon,...
      'Algorithm',                 'sqp',...
      'Display',                   'iter',...
      'FunctionTolerance',         1.0e-4,...
      'OptimalityTolerance',       1.0e-4);

    [qn_opt, fval] = fmincon(@(qq)obj_q(qq, S(:,n_t), K, ...
        call_P(:,n_t), put_P(:,n_t), B0(n_t)),...
        x0,[],[],Aeq,Beq,LB,UB,[],options);
    q(:, n_t) = qn_opt;
end

% Calculate transition probabilities p(i,j,n)
p0 = 1/m*ones(m*m,1);
p = zeros(m,m,N-1);
for n_t = 1:N-1
    Aeq_Metrix = [ones(1,m)];
    Aeq_log = repmat('Aeq_Metrix,',1,m);
    Aeq = eval(sprintf('blkdiag(%s)',Aeq_log(1:end-1)));
    Beq = [ones(1,m)];
    Beq = reshape(Beq, m,1);
    LB = zeros(m*m, 1);
    UB = ones(m*m, 1);
    
    [p_opt, fval] = fmincon(@(pp)obj_p_givenq(pp, q(:, n_t),q(:, n_t+1)),...
        p0, [], [], Aeq, Beq, LB, UB, [], options); 
    p_opt = (reshape(p_opt,m,m))';
    p(:,:,n_t) = p_opt;
end
end


