function [f_eff, f_eff_dot,f_int,f_env]= f_eff_fcn(Fs,f_int,L ,f_env, fiber, f_eff_dot, f_eff, Ueff,Af)

Ur = 0.8; % activation level at which all the motor units are recruited (Song et al. 2008)
US_th = 0.001; % threshold for slow-twitch fiber
UF_th = Ur*0.6; % threshold for fast-twitch fiber



if fiber== 'slow'
    T_f1=34.3;
    T_f2=22.7;
    T_f3=47;
    T_f4=25.2;
    % activation-frequency relationship (Brown and Loeb 2000
    f_half = 8.5; % frequency at which the motor unit produces half of its maximal isometric force 
    fmin = 0.5*f_half; % minimum firing frequency of slow-twitch fiber
    fmax = 2*f_half; % maximum firing frequency of slow-twitch fiber
    fenv = (fmax-fmin)/(1-US_th).*(Ueff-US_th)+fmin;
    f_env = fenv/f_half;
    
    if f_eff_dot>=0
        T_f=T_f1*L^2+ T_f2* f_env;
    else
        T_f=(T_f3+T_f4*Af)/L;
    end
    
elseif fiber== 'fast'
    T_f1=20.6;
    T_f2=13.6;
    T_f3=28.2;
    T_f4=15.1;
    
    f_half=34;
    fmin=0.5*f_half;
    fmax=2*f_half;
    fenv = (fmax-fmin)/(1-UF_th).*(Ueff-UF_th)+fmin;
    f_env = fenv/f_half;
    
    if f_eff_dot>=0
        T_f=T_f1*L^2+ T_f2* f_env;
    else
        T_f=(T_f3+T_f4*Af)/L;
    end
end




% intermediate firing frequency of second-order excitation dynamics
% of slow-twitch fiber (f_half)
fint_dot = (f_env - f_int)/T_f;
f_int = fint_dot*1/Fs + f_int;

if f_int <0
    f_int=0;
end
if f_env<0
    f_env=0;
end

% effective firing frequency of slow-twitch fiber (f_half)
f_eff_dot = (f_int - f_eff)/T_f;
f_eff = f_eff_dot*1/Fs + f_eff;
        
if f_eff < 0
    f_eff = 0;
elseif f_eff>1
    f_eff=1;
end

end