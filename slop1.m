%% 模拟数据
a = 0; b = 0; c = 1; tau = 0.5;
x = 0.01 : 0.01 : 1;% 不能从 0 开始, 向量化有些小问题.
xx = [x;x];
xx(1,xx(1,:)>tau) = 0;xx(2,xx(2,:)<=tau) = 0;
slope = [a, c];
tmp = xx;tmp(tmp > 0) = 1;
intercept = [b, a * tau + b - c * tau] * tmp;
y = slope * xx + intercept;
plot(x, y)
%% 估计