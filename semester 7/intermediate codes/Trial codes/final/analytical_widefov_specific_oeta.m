clc
clear all;
close all;

                                                              %%%% width of x axis
lowcap=110; 
incre=2.5;
upcap=135;
wid=((upcap-lowcap)/incre)+1;


gamma=3.1;                                      
BER=zeros(wid,2);
h=zeros(wid,2);
ebno=zeros(wid,2);
theta=zeros(wid,2);
d=2.5;
snr=zeros(wid,1);

theta(2)=pi/3;
theta(1)=(pi/6);


for gg=1:1:2
    h(gg)=((gamma+1)/2)*(9.55e-5)*cos(theta(gg))*((9+(d^2))^(-1*((gamma+1)/2)));
    check=0;
for s=lowcap:incre:upcap
    check=check+1;
    snr(check,1)=s;
    ebno(gg)= (h(gg)^2) * 10^(s/10);
    ggn=erfc(sqrt(ebno(gg)*0.5))*0.5;
    BER(check,gg)=ggn;
end  
end



semilogy(snr(:,1),BER(:,1),'--+b');
hold on
semilogy(snr(:,1),BER(:,2),'-.or');
hold on;
legend('\theta=30 degree','\theta=60 degree');
grid on