
clc;
clear all;
close all;

L=[0.5 1 1.3];
V=-5:0.01:5;

F_slow=zeros(length(L),length(V));
F_fast=zeros(length(L),length(V));


for i=1:length(L)
    for j=1:length(V)
        if V(j)<=0
            F_slow(i,j)=FVcon_Slow(L(i),V(j));
            F_fast(i,j)=FVcon_fast(L(i),V(j));
        elseif V(j)>0
            F_slow(i,j)=FVecc_slow(L(i),V(j));
            F_fast(i,j)=FVecc_fast(L(i),V(j));
        end
    end
end



plot(V,F_slow(1,:));

hold on
plot(V,F_fast(1,:));
