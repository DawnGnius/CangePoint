function [y] = generateData_Jump(n, gamma, tau, sigma)
%Input   n: Number of data points, gamma: Jump size, tau: Jump position
%Output  y: Data points
x = 1/n : 1/n : 1; 
Jump_Num = length(tau);
for ii = 1 : Jump_Num
    if ii == 1
        xx = x;
        xx(xx > tau(ii)) = 0;
    else
        xx = [xx; x];
        xx(ii,xx(ii,:) > tau(ii)) = 0;
        xx(ii,xx(ii,:) <= tau(ii - 1)) = 0;
    end
end % generate xx, accroding to tau
xx_01 = xx; xx_01(xx>0) = 1;
intercept = gamma * xx_01;
epsilon = sigma * randn(1,n);
y = 0 * xx + intercept + epsilon;
% plot(x,y,'.')