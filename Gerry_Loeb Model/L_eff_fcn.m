
function [Leff]= L_eff_fcn(Fs ,L, Af, L_eff0)

T_L= 0.088; %ms

L_eff_dot= (L-L_eff0)^3/T_L/(1-Af);
Leff= L_eff_dot*1/Fs +L_eff0;

end