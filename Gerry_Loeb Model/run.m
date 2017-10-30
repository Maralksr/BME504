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
L_eff_s=zeros(1,length(t)+1);
L_eff_f=zeros(1,length(t)+1);
S=zeros(1,length(t)+1);
f_int_slow=zeros(1,length(t)+1);
f_int_fast=zeros(1,length(t)+1);
f_env_fast=zeros(1,length(t)+1);
f_env_slow=zeros(1,length(t)+1);
f_eff_dot_slow=zeros(1,length(t)+1);
f_eff_dot_fast=zeros(1,length(t)+1);
f_eff_slow=zeros(1,length(t)+1);
f_eff_fast=zeros(1,length(t)+1);
Y=zeros(1,length(t)+1);
Af_slow=zeros(1,length(t)+1);
Af_fast=zeros(1,length(t)+1);
F_m=zeros(1,length(t)+1);
F_t=zeros(1,length(t)+1);

%L=[ones(1,Fs) 0.5*t(Fs+1:2*Fs)+0.5 1.5*ones(1,length(t)-2*Fs)]; 
%Lse=2-L;
%V=[zeros(1,Fs) 0.5*ones(1,Fs) zeros(1,length(t)-2*Fs)]; 

L=[ones(1,2*Fs) 1.5*ones(1,length(t)-2*Fs)]; 
Lse=2-L;
V=[zeros(1,2*Fs) 5 zeros(1,length(t)-2*Fs-1)]; 



Ueff_dot=zeros(1,length(t)+1);
Ueff=zeros(1,length(t)+1);
U=zeros(1,length(t)+1);
for i=1:length(t)
    
 %active force
 if U(i)>=Ueff(i)
     Tu=0.03;
 elseif U(i)<Ueff(i)
     TU=0.15;
 end
  
Ueff_dot(i)=(U(i)-Ueff(i))/Tu;
Ueff(i+1)=Ueff_dot(i)*1/Fs+U(i);
 
[Leff_s]= L_eff_fcn(Fs,L(i) ,Af_slow(i), L_eff_s(i));
L_eff_s(i+1)=Leff_s;

[Leff_f]= L_eff_fcn(Fs,L(i) ,Af_fast(i), L_eff_f(i));
L_eff_f(i+1)=Leff_f;
 
[S]=Sag_fcn(Fs, f_eff_slow(i), S(i));
S(i+1)=S;

[f_eff_slow, f_eff_dot_slow ,f_int_slow, f_env_slow]= f_eff_fcn(Fs, f_int_slow(i), L(i), f_env_slow(i), 'slow', f_eff_dot_slow(i), f_eff_slow(i), Ueff(i), Af_slow(i));
[f_eff_fast, f_eff_dot_fast ,f_int_fast, f_env_fast]= f_eff_fcn(Fs, f_int_fast(i), L(i), f_env_fast(i), 'fast', f_eff_dot_fast(i), f_eff_fast(i), Ueff(i), Af_fast(i));

f_int_slow(i+1)=f_int_slow;
f_int_fast(i+1)=f_int_fast;
f_eff_dot_slow(i+1)=f_eff_dot_slow;
f_eff_dot_fast(i+1)=f_eff_dot_fast;
f_eff_slow(i+1)=f_eff_slow;
f_eff_fast(i+1)=f_eff_fast;
f_env_slow(i+1)=f_env_slow;
f_env_fast(i+1)=f_env_fast;




[Y]= Yield_fnc(Fs, V(i), Y(i));
Y(i+1)=Y;

[Af_slow, Af_fast]=activ_force(Y(i),S(i),L_eff_s(i),L_eff_f(i), f_eff_slow(i), f_eff_fast(i));
Af_slow(i+1)=Af_slow;
Af_fast(i+1)=Af_fast;

%force length
FL_f=FL_fast(L(i));
FL_s=FL_slow(L(i));

if V(i)<=0
    FV_s=FVcon_slow(L(i),V(i));
    FV_f=FVcon_fast(L(i),V(i));
else
    FV_s=FVecc_slow(L(i),V(i));
    FV_f=FVecc_fast(L(i),V(i));
end

Ur = 0.8; % activation level at which all the motor units are recruited (Song et al. 2008)
US_th = 0.001; % threshold for slow-twitch fiber
dif_U_S = Ueff(i) - US_th; % difference between effective neural drive and threshold for slow-twitch fiber
UF_th = Ur*0.6; % threshold for fast-twitch fiber
dif_U_F = Ueff(i) - UF_th; % difference between effective neural drive and threshold for fast-twitch fiber
        if dif_U_F < 0
            dif_U_F = 0;
        end
        
        
W1 = dif_U_S/(dif_U_S+dif_U_F); % proportion of active slow-twitch fiber of total active muscle (0-1)
W2 = dif_U_F/(dif_U_S+dif_U_F); % proportion of active fast-twitch fiber of total active muscle (0-1)

FPE1=Fpe1(L(i),V(i)); %Stretching contractile passive element (fascicle) force (F0)
FPE2=Fpe2(L(i));

%----------
%total force

FCE_f = FL_f*FV_f + FPE2;
FCE_s = FL_s*FV_s + FPE2;
Fce= Ueff(i)*(W1.*Af_slow.*FCE_s+W2.*Af_fast.*FCE_f); %norm active force from contractile element with 2 fiber types
       
% total force from contractile element
F_muscle = Fce(1) + FPE1;


%tendon force
F_tendon=Fse(Lse(i));


F_m(i+1)=F_muscle;
F_t(i+1)=F_tendon;
MT_force(i+1)=F_muscle+ F_tendon;

i

end
