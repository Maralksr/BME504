function [Pulse]=Signal_Function(tp,ti,t) %This is the input signal which is defined as a function of t.
tpulse=ti+tp;   %Period ms
u=t./tpulse-floor(t./tpulse);   %Decimal
Pulse=(-sign(u-tp/tpulse)+1)/2;   %Using sign function to create 0 and 1
end