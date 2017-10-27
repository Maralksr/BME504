function [Lce_initial,Lse_initial,Lmax] =  InitialLength(Input)
        %---------------------------
        % Determine the initial lengths of muscle and tendon and maximal
        % muscle length 
        %---------------------------
        
        % serires elastic element parameters 
        cT = 27.8;
        kT = 0.0047;
        LrT = 0.964;
        % parallel passive element parameters 
        c1 = 23;
        k1 = 0.046;
        Lr1 = 1.17;
        
        % passive force produced by parallel passive element at maximal
        % muscle length
        PassiveForce = c1 * k1 * log(exp((1 - Lr1)/k1)+1);
        
        % tendon length at the above passive force 
        Normalized_SE_Length = kT*log(exp(PassiveForce/cT/kT)-1)+LrT;
        
        % maximal musculotendon length defined by joint range of motion
        Lmt_temp_max = Input.optimalLength*cos(Input.pennationAngle) ...
            +Input.tendonLength + 1;
        
        % optimal muscle length 
        L0_temp = Input.optimalLength;
        % optimal tendon length (Song et al. 2008)
        L0T_temp = Input.tendonLength*1.05;
        
        % tendon length at maximal muscle length
        SE_Length =  L0T_temp * Normalized_SE_Length;
        % maximal fasicle length
        FasclMax = (Lmt_temp_max - SE_Length)/L0_temp;
        % maximal muscle fiber length
        Lmax = FasclMax/cos(Input.pennationAngle);
        
        % initial musculotendon length defined by the user input
        Lmt_temp = Input.muscleInitialLength * cos(Input.pennationAngle) + Input.tendonInitialLength;
        
        % initial muscle length determined by passive muscle force and
        % tendon force
        InitialLength =  (Lmt_temp-(-L0T_temp*(kT/k1*Lr1-LrT-kT*log(c1/cT*k1/kT))))/(100*(1+kT/k1*L0T_temp/Lmax*1/L0_temp)*cos(Input.pennationAngle));
        % normalize the muscle legnth to optimal muscle length
        Lce_initial = InitialLength/(L0_temp/100);
        % calculate initial length of tendon and normalize it to optimal
        % tendon length
        Lse_initial = (Lmt_temp - InitialLength*cos(Input.pennationAngle)*100)/L0T_temp;
    end