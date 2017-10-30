clc;
clear all
close all
% initialization of activation dynamics parameters
Y_dot = 0;
S_dot_init = 0;
S_init = 0;
fint_dot = 0;



fint_dot_2 = 0;
fint_2 = 0;
feff_dot_2 = 0;
feff_2 = 0;


Fs=10000;
t=0:1/Fs:5;


Af=0;
L_init=0;




%%

%%
MT_force=zeros(1,length(t)+1);
L_eff=zeros(1,length(t)+1);
S=zeros(1,length(t)+1);
f_int_slow=zeros(1,length(t)+1);
f_int_fast=zeros(1,length(t)+1);
f_env_slow=zeros(1,length(t)+1);
f_env_fast=zeros(1,length(t)+1);
f_eff_dot_slow=zeros(1,length(t)+1);
f_eff_dot_fast=zeros(1,length(t)+1);
f_eff_slow=zeros(1,length(t)+1);
f_eff_fast=zeros(1,length(t)+1);
Y=zeros(1,length(t)+1);
Af_slow=zeros(1,length(t)+1);
Af_fast=zeros(1,length(t)+1);
F_m=zeros(1,length(t)+1);
F_t=zeros(1,length(t)+1);

L=1.5; 
Lse=1;
V=2;
Ueff=[0.2*t(1:length(t))];



for i=1:length(t)
    
 %active force
 [L_eff, L_eff_dot]=L_eff_fcn(Fs,L ,Af, L_eff(i));
 L_eff(i+1)=L_eff;
 
[S, S_dot]=Sag_fcn(Fs, f_eff_slow, S(i));
S(i+1)=S;

[f_eff_slow, f_eff_dot_slow ,f_int_slow, f_env_slow]= f_eff_fcn(Fs, f_int_slow(i), L, f_env_slow(i), 'slow', f_eff_dot_slow(i), f_eff_slow(i), Ueff(i), Af_slow(i));
[f_eff_fast, f_eff_dot_fast ,f_int_fast, f_env_fast]= f_eff_fcn(Fs, f_int_fast(i), L, f_env_fast(i), 'fast', f_eff_dot_fast(i), f_eff_fast(i), Ueff(i), Af_fast(i));
f_int_slow(i+1)=f_int_slow;
f_int_fast(i+1)=f_int_fast;
f_env_slow(i+1)=f_env_slow;
f_env_fast(i+1)=f_env_fast;
f_eff_dot_slow(i+1)=f_eff_dot_slow;
f_eff_dot_fast(i+1)=f_eff_dot_fast;
f_eff_slow(i+1)=f_eff_slow;
f_eff_fast(i+1)=f_eff_fast;




[Y, Y_dot]= Yield_fnc(Fs, V, Y(i));
Y(i+1)=Y;

[Af_slow, Af_fast]=activ_force(Y(i),S(i),L_eff(i),f_eff_slow(i), f_eff_fast(i));
Af_slow(i+1)=Af_slow;
Af_fast(i+1)=Af_fast;

%force length
FL_f=FL_fast(L);
FL_s=FL_Slow(L);

if V<=0
    FV_slow=FVcon_Slow(L,V);
    FV_fast=FVcon_fast(L,V);
else
    FV_slow=FVecc_slow(L,V);
    FV_fast=FVecc_fast(L,V);
end

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

FPE1=Fpe1(L,V);
FPE2=Fpe2(L);

%----------
%total force

FCE_f = FL_f*FV_fast + FPE2;
FCE_s = FL_s*FV_slow + FPE2;
Fce= Ueff(i)*(W1.*Af_slow.*FCE_s+W2.*Af_fast.*FCE_f);
       
% total force from contractile element
F_muscle = Fce(1) + FPE1;


%tendon force
F_tendon=Fse(Lse);


F_m(i+1)=F_muscle;
F_t(i+1)=F_tendon;
MT_force(i+1)=F_muscle+ F_tendon;

i

end
