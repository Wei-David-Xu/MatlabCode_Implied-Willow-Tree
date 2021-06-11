function call_Asian = Price_AsianCall_WT(T, r, S, p, q, S0, K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%   Calculate American call and put options prices by willow tree method.
%
% Input:
%   T : maturity
%   r : risk-free interest rate
%   S: willow tree of stock price -- m*N
%   p: transition probabilities -- m*m*N
%   q: transition probability from t0 to t1 -- m*1
%   S0: initial stock price
%   K: strike prices  -- n_K*1 
%
% Output:
%   call_Asian: vector of prices of a Asian call option of European style
%
% Implemented by
%      Bing Dong at April 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m, N] = size(S);
k_ave = 20;
Num_K = length(K);
% pricing Asian option for each strike K
for i = 1:Num_K
    [call_Asian(i)] = WTAsianE_oneStrike(m, N, T, r, S, p, q,S0, K(i),k_ave);
end
end

function [V] = WTAsianE_oneStrike(m, N, T, r, S, P, q, S0, K, k_ave)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Introduction:
%     This function returns the price of a Asian call option of European style 
%     by linear interpolation in Willow Tree method
% input:
%       m: number of spatial nodes in one moment
%       N: number of steps in Willow Tree, n+1 observed time points
%       T: the maturity of the Asian option
%       r: risk free interest rate
%       S: willow tree of stock price -- m*N
%       P: transition probabilities -- m*m*N
%       q: transition probability from t0 to t1 -- m*1
%       S0: the price of the underlying asset at time 0
%       K: the strike price
%       k_ave: number of average price in one moment
% output:
%       V: the price of a Asian call option of European style
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = T/N;
% max and min of average price in every moment
range = zeros(N, 2);
for i=1:N
    range(i,1) = (S0+sum(S(1,1:i)))/(i+1);
    range(i,2) = (S0+sum(S(m,1:i)))/(i+1);
end

% compute the option value at maturity
for k = 1:k_ave
    laststates(k) = range(N,1)+(k-1)*(range(N,2)-range(N,1))/(k_ave-1);
end

% compute the option value before maturity
 lastvalues = max(repmat(laststates, m, 1) - K, 0);
% lastvalues = max(K-repmat(laststates, m, 1) , 0);
 
for i = N-1:-1:1
    hh = (range(i + 1, 2)-range(i + 1, 1))/(k_ave-1);
    % present states and option values considered
    presentstates = zeros(1,k_ave);
    presentvalues = zeros(m,k_ave);
    for k = 1:k_ave
        presentstates(k) = range(i,1)+(k-1)*(range(i,2)-range(i,1))/(k_ave-1); 
    end
    for k=1:k_ave
        % avesum is used for interpolation at next time step
        avesum(:,k) = (presentstates(k)*(i+1)+S(:,i+1))/(i+2); 
        index = (avesum(:,k) - range(i + 1, 1)) / hh; 
        for kk = 1:m
            if index(kk) < 0
                index(kk) = 0;
            end
        end
        indexfloor = floor(index);
        interindex = indexfloor + (indexfloor~= k_ave-1);  
        % get xi, yi for interpolation
        x1 = laststates(interindex)';
        x2 = laststates(interindex + 1)';
        interindex = interindex + (0: k_ave : (m - 1) * k_ave)';
        lastvalues1 = lastvalues';
        y1 = lastvalues1(interindex);
        y2 = lastvalues1(interindex + 1);
        intervalue(:,k) = y1 + (avesum(:,k) - x1).*(y2 - y1)./ (x2 - x1);% linear interpolation
        presentvalues(:,k) = exp(-r*dt)*P(:,:,i)*intervalue(:,k);
    end
    laststates = presentstates;
    lastvalues = presentvalues;
end
% compute the option value at time 0
hh = (range(1, 2)-range(1, 1))/(k_ave-1);
avesum0 = (S0 + S(:, 1)) / 2;
index = (avesum0 - range(1, 1)) / hh;
indexfloor = floor(index);
interindex = indexfloor + (indexfloor~= k_ave-1);  
x1 = laststates(interindex)';
x2 = laststates(interindex + 1)';
interindex = interindex + (0: k_ave : (m - 1) * k_ave)';
lastvalues = lastvalues';
y1 = lastvalues(interindex);
y2 = lastvalues(interindex + 1);
intervalue0 = y1 + (y2 - y1) ./ (x2 - x1) .* (avesum0 - x1);% linear interpolation
V = exp(-r * dt) *q' * intervalue0;

end