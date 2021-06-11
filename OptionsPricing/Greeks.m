function [delta,gamma,theta] = Greeks(nodes,V,S0,V0,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Introdutcion
%         Given the willow tree including tree-nodes matrix and value
%          matrix, this function calculate the Greeks of an option,using the
%          Taylor expension.
%
%      Input
%         nodes (M*N matrix)         willow tree-nodes matrix
%          V  (M*N matrix)             V£¨i,j) represents the value of an option
%                                             at time j and node i
%           T(a scalar)                          expiration
%           S0 (a scalar)                       initial stock price
%           V0 (a scalar)                       the price of the option
%
%       Output
%          delta                   the delta of the option
%          gamma                   the gamma of the option
%          theta                   the theta of the option
%
%   Implememted by
%                P.Ai 2018.9.26
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,N] = size(nodes);
h = T/N;
e = ones(m,1);
A = [nodes(:,1)-S0 e*h ((nodes(:,1)-S0).^2)/2 (nodes(:,1)-S0)*h e*(h^2)/2];
b = (V(:,1)-V0*e);
x = A\b;
delta = x(1);
gamma = x(3);
theta = x(2);
end