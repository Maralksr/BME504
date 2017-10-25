
    function FV_ecc = FVecc_fast(L,V)
        %---------------------------
        % eccentric force velocity (F-V) relationship for fast-twitch fiber
        % input: normalized muscle length and velocity
        % output: F-V factor (0-1)
        %---------------------------
        av0 = -1.53;
        av1 = 0;
        av2 = 0;
        bv = 0.69;
        FV_ecc = (bv - (av0 + av1*L + av2*L^2)*V)/(bv+V);
    end
