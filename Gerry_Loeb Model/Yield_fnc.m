function [Y]= Yield_fnc(Fs, V, Y)
%for slow

C_Y=0.35;
V_Y=0.1;
T_Y= 200; %ms

Y_dot= (1-C_Y* (1-exp(-abs(V/V_Y)))-Y)/T_Y;
Y=Y_dot*1/Fs+Y;
end

