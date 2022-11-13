clc
clear all
peu=zeros(50,1);
pen=zeros(50,1);
snr=zeros(50,1);
l=4e-12;

for s=1:1:50
    k=s+109;
    en=10^(k/10);
    snr(s,1)=k;
    
    run=0;
    for x=l/2000:l/2000:l/2
        
        f1=(2/(3.147));                       
        f2=sqrt((4*x)*(l-x));
        f=f1/f2;
        ggu=erfc(sqrt(en*x/2))*0.5*f;             % for uniform distribution
        peu(s,1)=peu(s,1)+ ggu*(l/2000);

        run=run+1;
        
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(l-x));                  % for normal distribution
        ff3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        ff=2.85/(ff1*ff2*ff3);
        ggn=erfc(sqrt(en*x/2))*0.5*ff;
        pen(s,1)=pen(s,1)+ ggn*(l/2000);
        
    end
end
semilogy(snr,peu,'-+g');
hold on;
semilogy(snr,pen,'-+r');  
grid;
title('BER.VS.SNR)');
ylabel('BER');
xlabel('SNR [dB]');
legend('uniform','normal');