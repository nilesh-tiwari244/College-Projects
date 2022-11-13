clc;
clear all;
close all;

N=10000;
width=35;
lower_limit=80;
incre=2.5;
gamma=1;
x=randi([0 1],1,N);
mod_data=pskmod(x,2);
BER=zeros(2,width/incre);
new_data=zeros(2,N);
check1=0;

for SNR=lower_limit:incre:lower_limit+width
    
check1=check1+1;
err1=0;
err2=0;
jjjjj=0;
en=1/(10^(SNR/10));

for tt=1:1:N
    
 d=2;
 oeta1=pi/6;
 oeta2=pi/3;
        h1=((gamma+1)/2)*(9.55e-5)*cos(oeta1)*((9+(d^2))^(-1*((gamma+1)/2)));
        h2=((gamma+1)/2)*(9.55e-5)*cos(oeta2)*((9+(d^2))^(-1*((gamma+1)/2)));
noi=normrnd(0,sqrt(en),[1 N]);
recent=mod_data(tt);
new_data(1,:)= noi + h1*(mod_data);
new_data(2,:)= noi + h2*(mod_data);

if (new_data(1,tt)>0) && (recent<0)
    err1=err1+1;
elseif (new_data(1,tt)<0) && (recent>0)
       err1=err1+1;
   else 
       err1=err1+0;
   end


if (new_data(2,tt)>0) && (recent<0)
    err2=err2+1;
elseif (new_data(1,tt)<0) && (recent>0)
       err2=err2+1;
   else 
       err2=err2+0;
end

jjjjj=jjjjj+1;
check2=jjjjj+(check1)*i

end

BER(1,check1)=(err1/N);
BER(2,check1)=(err2/N);
end

SNR=lower_limit:incre:lower_limit+width;
ylabel('BER');
xlabel('SNR [dB]');
semilogy(SNR,BER(1,:),'-+gr');
hold on;
semilogy(SNR,BER(2,:),'-or');
legend('\theta=30 degree','\theta=60 degree');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55 

