clc
pe=zeros(50,1);
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
        gg=erfc(sqrt(en*x/2))*0.5*f;
        run=run+1;
        pe(s,1)=pe(s,1)+ gg*(l/2000);
    end
end
semilogy(snr,pe,'-+g');
hold on;
grid;
title('BER.VS.SNR)');
ylabel('BER');
xlabel('SNR [dB]');