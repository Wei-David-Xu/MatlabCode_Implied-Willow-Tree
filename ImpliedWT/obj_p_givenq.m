function SE = obj_p_givenq(p, q1, q2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Objective function when optimize transition probabilities from t to
%   t+dt in the willow tree given probabilities of each nodes,
%   i.e., Mean-Square-Error of |q1'* p - q2'|
%
% Input:
%   p : transition probability matrix
%   q1: probabilities of each nodes at t in willow tree
%   q2: probabilities of each nodes at t+dt in willow tree
%
% Output:
%   MSE: Mean-Square-Error of |V_mkt - V_mdl|
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = length(q1);
p = (reshape(p,m,m))';
SE = sum((q1'* p - q2').^2);
end
