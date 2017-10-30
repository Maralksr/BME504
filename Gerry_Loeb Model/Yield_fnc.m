function [Y, Y_dot]= Yield_fnc(Fs, V, Y_init)
%for slow

C_Y=0.35;
V_Y=0.1;
T_Y= 200; %ms

Y_dot= (1-C_Y* (1-exp(-abs(V/V_Y)))-Y_init)/T_Y;
Y=Y_dot*1/Fs+Y_init;
end

