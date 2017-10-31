function [x, xdot]= MyFun(f,t);
xdot=zeros(1,length(t));
x=zeros(1,length(t));
g=9.81;
m=1/g;
for i=1:length(t)
    xdot(i)=integral(f/m-g,0,t(i));
    x(i)=integral(xdot,0,t(i));
end
end



