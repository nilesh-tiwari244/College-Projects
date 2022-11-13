
clc;
clear all;
close all;

N=10000;
width=60;
lower_limit=90;
incre=2.5;


x=randi([0 1],1,N);
mod_data=pskmod(x,2);
BER=zeros(3,width/incre);
new_data=zeros(3,N);

gamma=3;
d=2.5;
fov1=pi/6;
fov2=pi/3;

check1=0;

for SNR=lower_limit:incre:lower_limit+width
    
check1=check1+1;
err1=0;
err2=0;
err3=0;
check2=0;
en=1/(10^(SNR/10));
 
for tt=1:1:N
         
theta=(pi/2)*rand(1);
if abs(theta) < fov1
    h1=((gamma+1)/2)*(9.55e-5)*abs(cos(theta))*((9+(d^2))^(-1*((gamma+1)/2)));
else
    h1=0;
end

if abs(theta) < fov2
    h2=((gamma+1)/2)*(9.55e-5)*abs(cos(theta))*((9+(d^2))^(-1*((gamma+1)/2)));
else
    h2=0;
end

 h3=((gamma+1)/2)*(9.55e-5)*abs(cos(theta))*((9+(d^2))^(-1*((gamma+1)/2)));

noi=normrnd(0,sqrt(en));
recent=mod_data(tt);
new_data(1,:)= noi + h1*(recent);
new_data(2,:)= noi + h2*(recent);
new_data(3,:)= noi + h3*(recent);

if (new_data(1,tt)>0) && (recent<0)
    err1=err1+1;
elseif (new_data(1,tt)<0) && (recent>0)
       err1=err1+1;
   else 
       err1=err1+0;
   end

if (new_data(2,tt)>0) && (recent<0)
    err2=err2+1;
elseif (new_data(2,tt)<0) && (recent>0)
       err2=err2+1;
   else 
       err2=err2+0;
end

if (new_data(3,tt)>0) && (recent<0)
    err3=err3+1;
elseif (new_data(3,tt)<0) && (recent>0)
       err3=err3+1;
   else 
       err3=err3+0;
end
check2=check2+1;
check3=check2+(check1)*i


end

BER(1,check1)=(err1/N);
BER(2,check1)=(err2/N);
BER(3,check1)=(err3/N);

end

SNR=lower_limit:incre:lower_limit+width;
ylabel('BER');
xlabel('SNR [dB]');
semilogy(SNR(:),BER(1,:),'-og');
hold on
semilogy(SNR(:),BER(2,:),'-+b');
hold on;
semilogy(SNR(:),BER(3,:),'-*r');
legend('FOV=35 degree','FOV=60 degree','WIDE FOV');
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

