 % (Song et al. 2008)
 Ueff=0;
 Y=0;
 Fs=10000;
 feff_dot=0;
 L=1;
 fint=0;
 feff=0;
 Ur = 0.8; % activation level at which all the motor units are recruited (Song et al. 2008)
 U1_th = 0.001; % threshold for slow-twitch fiber
  dif_U_1 = Ueff - U1_th; % difference between effective neural drive and threshold for slow-twitch fiber
  U2_th = Ur*0.6; % threshold for fast-twitch fiber
        dif_U_2 = Ueff - U2_th; % difference between effective neural drive and threshold for fast-twitch fiber
        if dif_U_2 < 0
            dif_U_2 = 0;
        end
        
        W1 = dif_U_1/(dif_U_1+dif_U_2); % proportion of active slow-twitch fiber of total active muscle (0-1)
        W2 = dif_U_2/(dif_U_1+dif_U_2); % proportion of active fast-twitch fiber of total active muscle (0-1)
        
        % activation-frequency relationship (Brown and Loeb 2000)
        f_half = 8.5; % frequency at which the motor unit produces half of its maximal isometric force 
        fmin = 0.5*f_half; % minimum firing frequency of slow-twitch fiber
        fmax = 2*f_half; % maximum firing frequency of slow-twitch fiber
        % constants for slow-twitch fiber
        af = 0.56;
        nf0 = 2.11;
        nf1 = 5;
        cy = 0.35;
        Vy = 0.1;
        Ty = 0.2;
        
        Tf1 = 0.0343;
        Tf2 = 0.0227;
        Tf3 = 0.047;
        Tf4 = 0.0252;
        
        % Y = yielding factor for slow-twitch fiber
        Y_dot = (1 - cy*(1-exp(-abs(V)/Vy))-Y)./Ty; 
        Y = Y_dot*1/Fs + Y;
        
        % firing frequency input to second-order excitation dynamics of
        % slow-twitch fiber
        fenv = (fmax-fmin)/(1-U1_th).*(Ueff-U1_th)+fmin;
        fenv = fenv/f_half;
        
        if feff_dot >= 0
            Tf = Tf1 * L^2 + Tf2 * fenv; % time constant for second-order excitation dynamics
        elseif feff_dot < 0
            Tf = (Tf3 + Tf4*Af_slow)/L;
        end
        
        % intermediate firing frequency of second-order excitation dynamics
        % of slow-twitch fiber (f_half)
        fint_dot = (fenv - fint)/Tf;
        fint = fint_dot*1/Fs + fint;
        % effective firing frequency of slow-twitch fiber (f_half)
        feff_dot = (fint - feff)/Tf;
        feff = feff_dot*1/Fs + feff;
        
        if feff < 0
            feff = 0;
        end
        
        % activation-frequency relationship for slow-twitch fiber
        nf = nf0 + nf1*(1/L-1);
        Af_slow = 1 - exp(-((Y*feff/(af*nf)).^nf));
        
        
        f_half_2 = 34;% frequency at which the motor unit produces half of its maximal isometric force 
        fmin_2 = 0.5*f_half_2; % minimum firing frequency of fast-twitch fiber
        fmax_2 = 2*f_half_2; % maximum firing frequency of fast-twitch fiber
        % constants for fast-twitch fiber
        af_2 = 0.56;
        nf0_2 = 2.1;
        nf1_2 = 3.3;
        as1_2 = 1.76;
        as2_2 = 0.96;
        Ts_2 = 0.043;
        
        Tf1_2 = 0.0206;
        Tf2_2 = 0.0136;
        Tf3_2 = 0.0282;
        Tf4_2 = 0.0151;
        
        % firing frequency input to second-order excitation dynamics of
        % fast-twitch fiber
        fenv_2 = (fmax_2-fmin_2)/(1-U2_th).*(Ueff-U2_th)+fmin_2;
        fenv_2 = fenv_2/f_half_2;
        
        if feff_dot_2 >= 0
            Tf_2 = Tf1_2 * L^2 + Tf2_2 * fenv_2; % time constant for second-order excitation dynamics
        elseif feff_dot_2 < 0
            Tf_2 = (Tf3_2 + Tf4_2*Af_fast)/L;
        end
        
        % Sagging factor (Saf) for fast-twitch fiber
        if feff_2 < 0.1
            as_2 = as1_2;
        elseif feff_2 >= 0.1
            as_2 = as2_2;
        end
        
        Saf_dot = (as_2 - Saf)/Ts_2;
        Saf = Saf_dot*1/Fs + Saf;
        
        % intermediate firing frequency of second-order excitation dynamics
        % of fast-twitch fiber (f_half)
        fint_dot_2 = (fenv_2 - fint_2)/Tf_2;
        fint_2 = fint_dot_2*1/Fs + fint_2;
        % effective firing frequency of fast-twitch fiber (f_half)
        feff_dot_2 = (fint_2 - feff_2)/Tf_2;
        feff_2 = feff_dot_2*1/Fs + feff_2;
        
        if feff_2 < 0
            feff_2 = 0;
        end
        
        % activation-frequency relationship for fast-twitch fiber
        nf_2 = nf0_2 + nf1_2*(1/L-1);
        Af_fast = 1 - exp(-((Saf*feff_2/(af_2*nf_2))^nf_2));