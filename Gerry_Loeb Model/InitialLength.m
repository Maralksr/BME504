function [Lce_initial,Lse_initial,Lmax] =  InitialLength(muscle_parameter)
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
        Lmt_temp_max = muscle_parameter.optimalLength*cos(muscle_parameter.pennationAngle) ...
            +muscle_parameter.tendonLength + 1;
        
        % optimal muscle length 
        L0_temp = muscle_parameter.optimalLength;
        % optimal tendon length (Song et al. 2008)
        L0T_temp = muscle_parameter.tendonLength*1.05;
        
        % tendon length at maximal muscle length
        SE_Length =  L0T_temp * Normalized_SE_Length;
        % maximal fasicle length
        FasclMax = (Lmt_temp_max - SE_Length)/L0_temp;
        % maximal muscle fiber length
        Lmax = FasclMax/cos(muscle_parameter.pennationAngle);
        
        % initial musculotendon length defined by the user input
        Lmt_temp = muscle_parameter.muscleInitialLength * cos(muscle_parameter.pennationAngle) + muscle_parameter.tendonInitialLength;
        
        % initial muscle length determined by passive muscle force and
        % tendon force
        InitialLength =  (Lmt_temp-(-L0T_temp*(kT/k1*Lr1-LrT-kT*log(c1/cT*k1/kT))))/(100*(1+kT/k1*L0T_temp/Lmax*1/L0_temp)*cos(muscle_parameter.pennationAngle));
        % normalize the muscle legnth to optimal muscle length
        Lce_initial = InitialLength/(L0_temp/100);
        % calculate initial length of tendon and normalize it to optimal
        % tendon length
        Lse_initial = (Lmt_temp - InitialLength*cos(muscle_parameter.pennationAngle)*100)/L0T_temp;
    end