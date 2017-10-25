    function FL = FL_fast(L)
        %---------------------------
        % force length (F-L) relationship for fast-twitch fiber
        % input: normalized muscle length
        % output: F-L factor (0-1)
        %---------------------------
        beta = 1.55;
        omega = 0.75;
        rho = 2.12;
        
        FL = exp(-abs((L^beta - 1)/omega)^rho);
    end