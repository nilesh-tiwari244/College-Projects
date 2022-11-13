
% clc
% clear all;
% close all;
                                                            
lowcap=90; 
incre=1;
upcap=113;
wid=((upcap-lowcap)/incre)+1;

                                      
BERana=zeros(wid,1);
h=zeros(1,1);
ebno=zeros(wid,1);
theta=zeros(wid,1);

snrana=zeros(wid,1);


   % h=    6.6753e-06;  % single antenna
    h= 7.585e-06;  % double antenna
    check=0;
for s=lowcap:incre:upcap
    check=check+1;
    snrana(check,1)=s;
    ebno(s)= (h^2) * 10^(s/10);
    ggn=erfc(sqrt(ebno(s)*0.5))*0.5;
    BERana(check,1)=ggn;
end  
% figure
% ylabel('BER');
% xlabel('SNR [dB]');
% semilogy(snrana(:,1),BERana(:,1),'-g');