function [L_eff, L_eff_dot]= L_eff_fcn(Fs ,L, Af, L_eff_init)

T_L= 0.088; %ms

L_eff_dot= (L-L_eff_init)^3/T_L/(1-Af);
L_eff= L_eff_dot*1/Fs +L_eff_init;

end