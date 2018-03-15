%% ������� MUller �����ı�д, ������ͬ

%% ����
a = 0; b = 1; c= -0.5; n = 100; L = 10; tau = [0.25, 0.5, 1]; sigma = 0.5; repeat =100; 
x = 1/n : 1/n : 1;
xx = [x;x;x];
for ii = 1 : length(tau) 
    xx(ii,xx(ii,:)>tau(ii)) = 0; 
    if ii ~= 1
        xx(ii,xx(ii,:) <= tau(ii-1)) = 0; 
    end
end

xx(xx > 0) = 1;% ����ֻ���ǽؾ���, ���, ������, xx����Ҫ.
slope = [0, 0, 0];
intercept = [a, b, c] * xx;

%% ��¼����ר��
Step = 1; Begin = 3; End = 50; RangOfL = Begin : Step: End;
recorder_L=zeros(length(RangOfL),2);
recorder_re = zeros(repeat,2); %��¼�ظ�repeat�εĽ��, �õ�ƽ��ֵ��Ϊ���յĹ���ֵ.(�ù���ֵ�Ƿ����)
for L = RangOfL
    tmp = (n - L) ;
    A = [1; 1/tmp];
    for ii= 1+1 : L
        A = [A, [1; ii/tmp]];
    end
    A=A';
    for num = 1 : repeat
        epsilon = sigma * randn(1,n);
        y = slope * xx + intercept + epsilon;
    %     plot(x, y, '.')
        for ii = 1 : L
            tmp = 0;
            for jj = 1 : floor(n/2) - L
                tmp = tmp + (y(jj+ii) - y(jj))^2;
            end
            if ii == 1
                Z = tmp;
            else
                Z = [Z, tmp];
            end
        end
        Z=Z'./(n-L);
        beta = [0.5, 0; 0, 1]* (A' * A)^-1 * A' * Z;
        recorder_re(num,:) = beta';
    end
    recorder_L((L-Begin)/Step+1,:)=mean(recorder_re); %û�з�������var(recorder_re)
end
figure()
hold on
plot(RangOfL, recorder_L(:,1)) % ����
plot(RangOfL, recorder_L(:,2),'--') % ��Ծ��С