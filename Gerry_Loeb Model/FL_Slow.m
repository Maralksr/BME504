function FL = FL_slow(L)
        %---------------------------
        % force length (F-L) relationship for slow-tiwtch fiber
        % input: normalized muscle length
        % output: F-L factor (0-1)
        %---------------------------
        beta = 2.3;
        omega = 1.12;
        rho = 1.62;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
 end