clc
clear all;
close all;

dim=60;                                                     %%%% width of x axis
snrlowerlimit=90;                                           %%%% x axis is ber vs snr
parts=10000;                                                 %%%% larger it is better the results
incre=2.5;
hac=snrlowerlimit:incre:snrlowerlimit+dim;

BER=zeros(length(hac),3);
snr=zeros(length(hac),1);                                         

fov1=pi*(35/180);
fov2=pi*(60/180); 

delta=[0.1318 9.85e-12];                                       %%% 1----FOV=35 degree
product1=0.058161757735338 ;                                   %%% 2----FOV=60 degree
cdf1=0.441287994957038;                                        %%% 3----wide FOV
product2= 1.849314141758325e-12;
cdf2= 0.187747628604906;

cdf=[cdf1 cdf2];
fov=[fov1 fov2 ];
product=[product1 product2 0];
 const=[1.5069 1.3012 1.2926];                    %% from different program

gamma=3;                                 
d=2.5;
hc2=(((gamma+1)/2)*(9.55e-5)*((9+(d^2))^(-1*((gamma+1)/2))))^2;     %%   6.7451e-13                                                                    %  hc2= (hc)^2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for aa=1:1:2 
    count=0;
for s=snrlowerlimit:incre:snrlowerlimit+dim
    count=count+1;
    k=s;
    s=count;
    ebno=10^(k/10);
    en=1/ebno;
    snr(s,1)=k;

init=0;
fina=hc2;
chang=(fina-init)/parts
  for x=init+chang:chang:fina-chang
         zz=0.5*acos((2*x*(1/hc2))-1);
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc2-x));                  % for normal distribution with narrow fov
        ff3=exp(((zz-(pi/6))^2)/(2*(pi/9)));
        ff=(const(aa))/(ff1*ff2*ff3);
        ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        BER(s,aa)=BER(s,aa)+ ggn*chang;
    end
        x=0; 
        ffff=cdf(aa);
        ggggn=erfc(sqrt(ebno*x*0.5))*0.5*ffff;
       BER(s,aa)=BER(s,aa)+ ggggn;    
end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa=3;
count=0;
for s=snrlowerlimit:incre:snrlowerlimit+dim
    count=count+1;
k=s;
s=count;
ebno=10^(k/10);
snr(s,1)=k;
init=0;
fina=hc2;
chang=(fina-init)/parts;
  for x=init+chang:chang:fina-chang
         zz=0.5*acos((2*x*(1/hc2))-1);
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc2-x));                  % for normal distribution with wide FOV
        ff3=exp(((zz-(pi/6))^2)/(2*(pi/9)));
        ff=(const(aa))/(ff1*ff2*ff3);
         ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        BER(s,aa)=BER(s,aa)+ ggn*chang;
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

semilogy(snr,BER(:,1),'-g');
hold on
semilogy(snr,BER(:,2),'--b');
hold on
semilogy(snr,BER(:,3),'-.r');
grid on

ylabel('BER');
xlabel('SNR [dB]');
legend('FOV=35 degree','FOV=60 degree','WIDE FOV');



