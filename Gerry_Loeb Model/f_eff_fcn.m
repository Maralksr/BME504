function [f_eff, f_eff_dot,f_int, f_env]= f_eff_fcn(Fs,f_int,L,f_env, fiber, f_eff_dot_init, f_eff_init, Ueff_init,Af)

Ur = 0.8; % activation level at which all the motor units are recruited (Song et al. 2008)
U1_th = 0.001; % threshold for slow-twitch fiber
U2_th = Ur*0.6; % threshold for fast-twitch fiber



if fiber== 'slow'
    T_f1=24.2;
    T_f2=16;
    T_f3=33.2;
    T_f4=17.8;
    % activation-frequency relationship (Brown and Loeb 2000
    f_half = 8.5; % frequency at which the motor unit produces half of its maximal isometric force 
    fmin = 0.5*f_half; % minimum firing frequency of slow-twitch fiber
    fmax = 2*f_half; % maximum firing frequency of slow-twitch fiber
    fenv = (fmax-fmin)/(1-U1_th).*(Ueff_init-U1_th)+fmin;
    f_env = fenv/f_half;
    if f_eff_dot_init>=0
        T_f=T_f1*L^2+ T_f2* f_env;
    else
        T_f=(T_f3+T_f4*Af)/L;
    end
elseif fiber== 'fast'
    T_f1=20.6;
    T_f2=13.6;
    T_f3=28.2;
    f_half=34;
    fmin=0.5*f_half;
    fmax=2*f_half;
    fenv = (fmax-fmin)/(1-U2_th).*(Ueff_init-U2_th)+fmin;
    f_env = fenv/f_half;
    
    T_f=T_f1*L^2+ T_f2* f_env;
end








% intermediate firing frequency of second-order excitation dynamics
% of slow-twitch fiber (f_half)
fint_dot = (f_env - f_int)/T_f;
f_int = fint_dot*1/Fs + f_int;
% effective firing frequency of slow-twitch fiber (f_half)
f_eff_dot = (f_int - f_eff_init)/T_f;
f_eff = f_eff_dot*1/Fs + f_eff_init;
        
if f_eff < 0
    f_eff = 0;
end

end