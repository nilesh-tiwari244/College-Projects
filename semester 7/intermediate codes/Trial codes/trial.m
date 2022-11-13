clc
clear all;
close all;


dim=60;                                                     %%%% width of x axis
snrlowerlimit=80;                                           %%%% x axis is ber vs snr
parts=1000;                                                 %%%% larger it is better the results
                                         
BER=zeros(dim,3);
snr=zeros(dim,1);                                         

fov1=pi*(35/180);
fov2=pi*(60/180);                                              %%% 1----FOV=35 degree
                                                               %%% 2----FOV=60 degree
cdf1=0.6111;                                                   %%% 3----wide FOV
cdf2= 0.27777;
delta=[0.1318 9.85e-12];
product1=cdf1*delta(1) ;  
product2=cdf2*delta(2) ;  

cdf=[cdf1 cdf2];
fov=[fov1 fov2 ];
product=[product1 product2 0];
const=[3.7569 2.1270 1.4259];

gamma=1.2;                                 % gamma
d=2.5;
hc2=(((gamma+1)/2)*(9.55e-5)*((9+(d^2))^(-1*((gamma+1)/2))))^2         
%  hc2=4e-12                                                               %  hc2= (hc)^2                                                                 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for aa=1:1:2    
for s=1:1:dim
    k=s+snrlowerlimit;
    
    ebno=10^(k/10);
    en=1/ebno;
    snr(s,1)=k;
    for x=(hc2/parts)*(cos(fov(aa))^2):hc2/parts:hc2-hc2/parts
        ff2=sqrt((4*x)*(hc2-x));                                       % for uniform distribution with narrow fov
        ff=(2/pi)/(ff2);
        ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        BER(s,aa)=BER(s,aa)+ ggn*(hc2/parts);
    end
         x=0; 
        ffff=delta(aa);
        ggggn=erfc(sqrt(ebno*x*0.5))*0.5*ffff;
       BER(s,aa)=BER(s,aa)+ ggggn;    
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa=3;
for s=1:1:dim
    k=s+snrlowerlimit;
    ebno=10^(k/10);
    snr(s,1)=k;
    for x=hc2/parts:hc2/parts:hc2-hc2/parts
        ff2=sqrt((4*x)*(hc2-x));                  % for uniform distribution with wide fov
        ff=(2/pi)/(ff2);
        ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        BER(s,aa)=BER(s,aa)+ ggn*(hc2/parts);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

semilogy(snr,BER(:,1),'-+b');
hold on
semilogy(snr,BER(:,2),'-*g');
hold on
semilogy(snr,BER(:,3),'-or');
grid on

ylabel('BER');
xlabel('SNR [dB]');
legend('FOV=35 degree','FOV=60 degree','WIDE FOV');


