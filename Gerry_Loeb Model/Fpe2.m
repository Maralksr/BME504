
    function Fpe2 = Fpe2(L)
        %---------------------------
        % passive element 2 or thick element compression
        % input: normalized muscle length
        % output: passive element force (0-1)
        %---------------------------
        c2_pe2 = -0.02; 
        k2_pe2 = -21;
        Lr2_pe2 = 0.70; 
        
        Fpe2 = c2_pe2*exp((k2_pe2*(L-Lr2_pe2))-1);
        
    end