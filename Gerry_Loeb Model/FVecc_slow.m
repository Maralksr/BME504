    function FV_ecc = FVecc_slow(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for slow-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -4.7;
        av1 = 8.41;
        av2 = -5.34;
        bv = 0.35;
        FV_ecc = -(bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end