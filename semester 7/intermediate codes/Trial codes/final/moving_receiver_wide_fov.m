clc
clear all
close all
dim=20;
peu=zeros(dim,1);
pen=zeros(dim,1);
snr=zeros(dim,1);
l=4e-12;

for s=1:1:dim
    k=s+120;
    en=10^(k/10);
    snr(s,1)=k;
    
    run=0;
%     for x=l*0.586:l/2000:l*0.883
       for x=l/20000:l/20000:l-(l/20000) 
        f1=(2/(3.147));                       
        f2=sqrt((4*x)*(l-x));
        f=f1/f2;
        ggu=erfc(sqrt(en*x/2))*0.5*f;            % for uniform distribution
        peu(s,1)=peu(s,1)+ ggu*(l/20000);

    end
%       for x=l/2000:l/2000:l-(l/2000)   
%           
%         ff1=sqrt(2*pi*(pi/9));
%         ff2=sqrt((4*x)*(l-x));                  % for normal distribution
%         ff3=exp(((x-(pi/6))^2)/(2*(pi/9)));
%         ff=1.42/(ff1*ff2*ff3);
%         ggn=erfc(sqrt(en*x/2))*0.5*ff;
%         pen(s,1)=pen(s,1)+ ggn*(l/2000);
%         
%     end
end
semilogy(snr,peu,'-*g');
hold on;
semilogy(snr,pen,'-or');  
grid;
title('BER.VS.SNR)');
ylabel('BER');
xlabel('SNR [dB]');
legend('uniform distribution','normal distribution');