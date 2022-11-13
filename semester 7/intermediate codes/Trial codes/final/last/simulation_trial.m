clc;
clear all;
close all;

N=1000;
x=randi([0 1],1,N);
mod_data=pskmod(x,2);


width=25;
lower_limit=110;
incre=5;
gamma=3.1;


BER=zeros(3,width/incre);
new_data=zeros(3,N);

d=2.5;
theta=[0 pi/6 pi/3];
err=[0 0 0];
h=[0 0 0];

for mm=3:1:3
 check1=0;   
 h(mm)=((gamma+1)/2)*(9.55e-5)*cos(theta(mm))*((9+(d^2))^(-1*((gamma+1)/2)));
 
 

for SNR=lower_limit:incre:lower_limit+width    
    
check1=check1+1;
check2=0;
en=1/(10^(SNR/10));



for tt=1:1:N
    
noi=normrnd(0,sqrt(en/2),[1 N]);
new_data(mm,:)= noi + h(mm)*(mod_data);
recent=mod_data(tt);

if (new_data(mm,tt)>0) && (recent<0)
    err(mm)=err(mm)+1;
elseif (new_data(mm,tt)<0) && (recent>0)
       err(mm)=err(mm)+1;
   else 
       err(mm)=err(mm)+0;
end
check2=check2+1;
check3=mm+(check1)*i

end


BER(mm,check1)=(err(mm)/N);
end

end

SNR=lower_limit:incre:lower_limit+width;
ylabel('BER');
xlabel('SNR [dB]');
% semilogy(SNR,BER(1,:),'--+gr');
% hold on;
% semilogy(SNR,BER(2,:),'-.or');
% hold on
semilogy(SNR,BER(3,:),':ob');
legend('\theta=0 degree','\theta=30 degree','\theta=60 degree');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55 

