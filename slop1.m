format compact
%% 模拟数据
a = 0; b = 0; c = 3; n = 500; L = 50; tau = 0; sigma = 0.03; repeat =1000; 
x = 1/n : 1/n : 1;% 不能从 0 开始, 向量化有些小问题.
xx = [x;x];
xx(1,xx(1,:)>tau) = 0; xx(2,xx(2,:)<=tau) = 0;
slope = [a, c];
tmp = xx;tmp(tmp > 0) = 1;
intercept = [b, a * tau + b - c * tau] * tmp;
recorder = zeros(repeat,2);
for num = 1: repeat
    epsilon = sigma * randn(1,n);
    y = slope * xx + intercept + epsilon;
    % plot(x, y)
    %% 估计
    tmp = (floor(n/2)-L) * n^2;
    A = [1; 1/tmp];
    for ii= 1+1 : L
        A = [A, [1; ii/tmp]];
    end%数字太小了, 会不会有很大的影响
    for ii = 1 : L
        tmp = 0;
        for jj = 1 : floor(n/2) - L
            tmp = tmp + (y(2*jj+2*ii) - y(2*jj+2*ii-1) - y(2*jj) + y(2*jj-1))^2;
        end
        if ii == 1
            Z = tmp;
        else
            Z = [Z, tmp];
        end
    end
    A = A'; Z = Z'./ (floor(n/2)-L);
    beta = (A' * A)^-1 * A' * Z;
    recorder(num,:) = beta';
end

%% Analysis
mean(recorder)
var(recorder)