function [S]= Sag_fcn(Fs, f_eff, S_init)
%sag fast
a_s1=1.76;
a_s2=0.96;
Ts= 43; %ms

if f_eff <0.1
    a_s=a_s1;
else
    a_s=a_s2;
end
S_dot= (a_s- S_init)/Ts;
S= S_dot*1/Fs+S_init;
end

