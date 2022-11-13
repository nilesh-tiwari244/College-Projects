tic

clc;
close all;
% global N mod_data x width incre lower_limit no_receiver BER gamma d fov1 fov diff theta_mean theta_var ;


N=1000;
x=randi([0 1],1,N);
mod_data=pskmod(x,2);

width=50;
incre=2.5;
lower_limit=100;

no_receiver=20;

BER=zeros(3,width/incre);

gamma=3;
d=2.5;

fov1=pi*(60/180);
fov=fov1*ones(no_receiver,1);
diff=[0 ; (rand([no_receiver-1 1])-0.5)*20*(pi/180)];
% diff=[0; 5; 10; -5; -10]*(pi/180);

theta_mean=pi*(60/180);
theta_var=0.01;

check1=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for SNR=lower_limit:incre:lower_limit+width
    new_data=zeros(3,N);
   
    check1=check1+1

err1=0;
err2=0;
err3=0;

en=1/(10^(SNR/10));
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for tt=1:1:N
    
    
    kk=zeros(no_receiver,1);
    ww=ones(no_receiver,1);
    hh=zeros(no_receiver,1);   
    
theta=abs(normrnd(theta_mean,sqrt(theta_var)));

 for re=1:1:no_receiver
    if abs(theta-diff(re,1)) < fov(re,1)
        hh(re)=((gamma+1)/2)*(9.55e-5)*(cos(theta-diff(re,1)))*((9+(d^2))^(-1*((gamma+1)/2)));
    else
        hh(re)=0;
    end
end

  h3=((gamma+1)/2)*(9.55e-5)*(cos(theta))*((9+(d^2))^(-1*((gamma+1)/2)));
 
 noi=normrnd(0,sqrt(en));
 recent=mod_data(tt);
 new_data(1,tt)= noi + hh(1)*(recent);
 new_data(3,tt)= noi + h3*(recent);

if (new_data(1,tt)>0) && (recent<0)
    err1=err1+1;
elseif (new_data(1,tt)<0) && (recent>0)
       err1=err1+1;
   else 
       err1=err1+0;
end

   if (new_data(3,tt)>0) && (recent<0)
    err3=err3+1;
elseif (new_data(3,tt)<0) && (recent>0)
       err3=err3+1;
   else 
       err3=err3+0;
   end


sq_sum=sum(hh.^2);


for re=1:1:no_receiver
    noi=normrnd(0,sqrt(en));
    kk(re,1)= noi + hh(re,1)*(recent);
end


if (sq_sum==0)
    for re=1:1:no_receiver
        ww(re,1)=1/sqrt(no_receiver); 
    end   
else
    for re=1:1:no_receiver
        ww(re,1)=hh(re,1)/sqrt(sq_sum); 
    end
end


for re=1:1:no_receiver
    new_data(2,tt)=new_data(2,tt) + kk(re,1)*ww(re,1);
end


if (new_data(2,tt)>0) && (recent<0)
    err2=err2+1;
elseif (new_data(2,tt)<0) && (recent>0)
       err2=err2+1;
   else 
       err2=err2+0;
end

end

BER(1,check1)=(err1/N);
BER(2,check1)=(err2/N);
BER(3,check1)=(err3/N);

end
figure
SNR=lower_limit:incre:lower_limit+width;
ylabel('BER');
xlabel('SNR [dB]');
semilogy(SNR(:),BER(1,:),'-og');
hold on
semilogy(SNR(:),BER(2,:),'-+b');
hold on;
semilogy(SNR(:),BER(3,:),'-*r');
legend('FOV=60 degree','multiple FOV=60 degree','WIDE FOV');
grid on

errors=[err1 err2 err3]
sq_sum
hh
ww
last_theta=theta*(180/pi)

satu_prob=func_norm_prob(theta_mean,theta_var,fov1,inf)

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
