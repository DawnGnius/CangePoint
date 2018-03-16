format compact
%% Simulation
a = 0; b = 0; c = 3; n = 200; tau = 0.5; sigma = 0.2; repeat =1000; 
x = 1/n : 1/n : 1;% Start from 0 could lead to a problem in the following vectorization process.
xx = [x;x];
xx(1,xx(1,:)>tau) = 0; xx(2,xx(2,:)<=tau) = 0;
slope = [a, c];
tmp = xx;tmp(tmp > 0) = 1;
intercept = [b, a * tau + b - c * tau] * tmp;

% Prepare for recorder_L, recorder_re
Begin = 2; Step = 1; End = 50; RangOfL = Begin : Step : End;
recorder_L=zeros(length(RangOfL),2);
recorder_L_var=zeros(length(RangOfL),2);
recorder_re = zeros(repeat,2);
for L = RangOfL
    recorder_re = zeros(repeat,2);
    for re_num = 1 : repeat
        epsilon = sigma * randn(1,n);
        y = slope * xx + intercept + epsilon;
        % plot(x, y)
        %% Estimation
        tmp = (floor(n/2)-L) * n^2;
        A = [1; 1/tmp];
        for ii= 1+1 : L
            A = [A, [1; ii/tmp]];
        end % Caution: Numeric error
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
        A = A'; Z = Z' ./ (floor(n/2)-L);
        beta = (A' * A)^-1 * A' * Z;
        recorder_re(re_num,:) = beta';
    end
    recorder_L((L-Begin)/Step+1,:)=mean(recorder_re); % Ignore the variance, var(recorder_re), dim: 1 * 2
    recorder_L_var((L-Begin)/Step+1,:)=var(recorder_re); % Strange setting, Caution; ignore var
end
%% Plot
figure()
hold on
plot(RangOfL, recorder_L(:,1)) % beta(1) is sigma^2
plot(RangOfL, recorder_L(:,2),'--') % beta(2) is gamma
figure()
hold on
plot(RangOfL, recorder_L_var(:,1)) % beta(1) is sigma^2 -- varicance
plot(RangOfL, recorder_L_var(:,2),'--') % beta(2) is gamma -- jump size