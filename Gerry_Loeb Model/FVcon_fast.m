    function FV_con = FVcon_fast(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -9.15;
        cv0 = -5.7;
        cv1 = 9.18;
        
        FV_con = (Vmax-V)/(Vmax + (cv0 + cv1*L)*V);
    end