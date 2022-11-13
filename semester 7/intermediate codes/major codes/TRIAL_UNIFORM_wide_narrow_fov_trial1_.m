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
fov2=pi*(60/180);                                              %%% 1----FOV=35 degree
                                                               %%% 2----FOV=60 degree
cdf1=0.611;                                                   %%% 3----wide FOV
cdf2= 0.33;

cdf=[cdf1 cdf2];
fov=[fov1 fov2 ];
const=[0.1334 0.3042 0.4057];

gamma=3;                                 
d=2.5;
hc2=(((gamma+1)/2)*(9.55e-5)*((9+(d^2))^(-1*((gamma+1)/2))))^2                                                                      %  hc2= (hc)^2                                                                 



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
chang=(fina-init)/parts;
  for x=init+chang:chang:fina-chang
        ff2=sqrt((4*x)*(hc2-x));                                       % for uniform distribution with narrow fov
        ff=(2/pi)/(ff2);
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
    for x=hc2/parts:hc2/parts:hc2-hc2/parts
        ff2=sqrt((4*x)*(hc2-x));                  % for uniform distribution with wide fov
        ff=(2/pi)/(ff2);
        ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        BER(s,aa)=BER(s,aa)+ ggn*(hc2/parts);
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


