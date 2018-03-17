function estimate_Jump_1(y)
    n = length(y);
    Step = 1; Begin = 2; End = 50; RangOfL = Begin : Step: End;
    recorder_L=zeros(length(RangOfL),2);
    for L = RangOfL
        tmp = (n - L) ;
        A = [1; 1/tmp];
        for ii= 1 + 1 : L
            A = [A, [1; ii/tmp]];
        end
        A=A';
        % plot(x, y, '.')
        for ii = 1 : L
            tmp = 0;
            for jj = 1 : n - L
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
        recorder_L((L-Begin)/Step+1,:)=beta; % Strange setting, Caution; ignore var
    end
    figure()
    hold on
    plot(RangOfL, recorder_L(:,1)) % beta(1) is sigma^2 -- varicance
    plot(RangOfL, recorder_L(:,2),'--') % beta(2) is gamma -- jump size
end % Function