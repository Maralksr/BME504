d=fdesign.lowpass('N,Fc',4,30,1000);
designmethods(d)
Hd=design(d);
%fvtool(Hd)
output=filter(Hd, EMG);


Activation=output/ (max(output));
plot(Activation)
max(Activation)