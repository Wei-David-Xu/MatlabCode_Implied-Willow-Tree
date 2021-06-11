function [itype,ifault,Willow,H] = WTnodes_from_JohnsonCurve(G,N,M,gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   This function returns the random variable Willow 
%   based on the Johnson curve inverse transformation
% 
% Input:
%       G: fourth-order moment of asset return ln(St/S0)
%       M: number of nodes at each time in the willow tree
%       N: number of time steps in Willow Tree
%       gamma: parameter of willow tree (default: 0.6)
% Output:
%       itype: determine the form of function g
%       ifault
%       Willow: the random variable matrix based on the Johnson curve inverse transformation
%       H:      contains values of parameters a,b,c,d
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% compute z,q
for k = 1:M/2+1
   q(k) = (k-0.5)^gamma/M;
   q(M+1-k) = q(k);
end
qsum = sum(q);
q = q./qsum;
z(1) = norminv(q(1)/2,0,1);
for k = 2:M
   tmp = sum(q(1:k-1))+q(k)/2;
   z(k) = norminv(tmp,0,1);
end

% call f_hhh with C++ to determine the form of g and compute a,b,c,d
for i = 1:N
    mu = G(1,i);
    sd = sqrt(G(2,i));
    ka3 = G(3,i);
    ka4 = G(4,i);
    [a,b,d,c,itype(i),ifault(i)] = f_hhh(mu,sd,ka3,ka4);
% a standard normally distributed random variable can be converted into a random variable 
    if itype(i) == 1
        u = (z(:)-a)./b;
        gi = exp(u(:));
        x = c+d.*gi(:);
    elseif itype(i) == 2
        u = (z(:)-a)./b;
        gi = ((exp(u(:))-exp(-u(:)))./2);
        x = c+d.*gi(:);
    elseif itype(i) == 3
        u = (z(:)-a)./b;
        gi = 1./(1+exp(-u(:)));
        x = c+d.*gi(:);
    elseif itype(i) == 4
        u = (z(:)-a)./b;
        gi = u(:);
        x = c+d.*gi(:);
    end
    H(1,i) = a;
    H(2,i) = b;
    H(3,i) = c;
    H(4,i) = d;
    Willow(:,i) = x;
end
end

