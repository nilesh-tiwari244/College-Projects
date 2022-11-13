clc
clear all;
close all;

d=2.5;
bounds=6
hc2=zeros(length(-bounds:0.01:bounds),1);
k=0;
for gamma=-bounds:0.01:bounds
    k=k+1;
    hc2(k)=(((gamma+1)/2)*(9.55e-5)*((9+(d^2))^(-1*((gamma+1)/2))))^2;
end
gamma=-bounds:0.01:bounds;
semilogy(gamma,hc2,'--or')
grid on;