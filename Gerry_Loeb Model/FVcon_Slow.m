   function FV_con = FVcon_slow(L,V)
        %---------------------------
        % concentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        Vmax = -7.88;
        cv0 = 5.88;
        cv1 = 0;
        
        FV_con = (Vmax-V)/(Vmax + (cv0 + cv1*L)*V);
    end