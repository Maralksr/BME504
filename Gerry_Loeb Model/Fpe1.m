    function Fpe1 = Fpe1(L,V)
        %---------------------------
        % passive element 1 
        % input: normalized muscle length
        % output: passive element force (0-1)
        %---------------------------
        c1_pe1 = 23;
        k1_pe1 = 0.046;
        Lr1_pe1 = 1.17;
        eta = 0.01;
        
        Fpe1 = c1_pe1 * k1_pe1 * log(exp((L - Lr1_pe1)/k1_pe1)+1) + eta*V;
        
    end
