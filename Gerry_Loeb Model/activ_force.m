function [Af_slow, Af_fast]=activ_force(Y,S,L_eff,f_eff_slow,f_eff_fast)
    


a_f_s=0.56;
n_f0_s=2.1;
n_f1_s=5;
    
a_f_f=0.56;
n_f0_f=2.1;
n_f1_f=3.3;
    

n_f_s= n_f0_s + n_f1_s*(1/L_eff-1);
n_f_f= n_f0_f + n_f1_f*(1/L_eff-1);

Af_slow= 1- exp(-(Y*f_eff_slow/a_f_s/n_f_s)^n_f_s);
Af_fast= 1- exp(-(S*f_eff_fast/a_f_f/n_f_f)^n_f_f);
end
