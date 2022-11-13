
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%                                  SIMULATION WIDE FOV(GAMMA=3.18)


clc;
clear all;
close all;

N=10000;
width=25;
lower_limit=110;
incre=2.5;


x=randi([0 1],1,N);
mod_data=pskmod(x,2);
BER=zeros(2,width/incre);
new_data=zeros(2,N);

gamma=3.18;
d=2.5;
theta1=pi/6;
theta2=pi/3;


check1=0;

for SNR=lower_limit:incre:lower_limit+width
    
check1=check1+1;
err1=0;
err2=0;
jjjjj=0;
en=1/(10^(SNR/10));

for tt=1:1:N
    
        h1=((gamma+1)/2)*(9.55e-5)*cos(theta1)*((9+(d^2))^(-1*((gamma+1)/2)));
        h2=((gamma+1)/2)*(9.55e-5)*cos(theta2)*((9+(d^2))^(-1*((gamma+1)/2)));
        
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
semilogy(SNR(:),BER(1,:),'--+b');
hold on
semilogy(SNR(:),BER(2,:),'-.or');
hold on;
legend('\theta=30 degree','\theta=60 degree');
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                              ANALYTICAL WIDE FOV (GAMMA=3.1)


clc
clear all;
close all;
                                                             %%%% width of x axis
lowcap=110; 
incre=2.5;
upcap=135;
wid=((upcap-lowcap)/incre)+1;

                                      
BER=zeros(wid,2);
h=zeros(wid,2);
ebno=zeros(wid,2);
theta=zeros(wid,2);

snr=zeros(wid,1);

gamma=3.1;
d=2.5;
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

ylabel('BER');
xlabel('SNR [dB]');
semilogy(snr(:,1),BER(:,1),'--+b');
hold on
semilogy(snr(:,1),BER(:,2),'-.or');
hold on;
legend('\theta=30 degree','\theta=60 degree');
grid on
 