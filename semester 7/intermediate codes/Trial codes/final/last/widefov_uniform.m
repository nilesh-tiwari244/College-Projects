clc
clear all;
close all;
dim=60;                                                     %%%% width of x axis
snrlowerlimit=80;                                             %%%% x axis is ber vs snr
parts=1000;                                                 %%%% larger it is better the results
                                         
pen=zeros(dim,3);
snr=zeros(dim,1);                                         

fov1=pi*(35/180);
fov2=pi*(60/180);                                              %%% 1----FOV=35 degree
                                 %%% 2----FOV=60 degree
cdf1=0.6111;                                        %%% 3----wide FOV
cdf2= 0.27777;
delta=[0.1318 9.85e-12];
product1=cdf1*delta(1) ;  
product2=cdf2*delta(2) ;  


cdf=[cdf1 cdf2];
fov=[fov1 fov2 ];
product=[product1 product2 0];
const=[3.7569 2.1270 1.4259];
hc3=4e-12;
hc=[hc3 hc3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for aa=1:1:2    
for s=1:1:dim
    k=s+snrlowerlimit;
    
    ebno=10^(k/10);
    en=1/ebno;
    snr(s,1)=k;
    for x=(hc(aa)/parts)*(cos(fov(aa))^2):hc(aa)/parts:hc(aa)-hc(aa)/parts
        % zz=0.5*acos((2*x*(1/hc(aa)))-1);
      %  ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc(aa)-x));                  % for normal distribution with narrow fov
      %  ff3=exp(((zz-(pi/6))^2)/(2*(pi/9)));
        ff=(2/pi)/(ff2);
        ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        %ggn=erfc(sqrt(ebno*x))*0.5*ff;
        pen(s,aa)=pen(s,aa)+ ggn*(hc(aa)/parts);
    end
         x=0; 
        ffff=delta(aa);
        ggggn=erfc(sqrt(ebno*x*0.5))*0.5*ffff;
       pen(s,aa)=pen(s,aa)+ ggggn;    
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aaa=3;
for s=1:1:dim
    k=s+snrlowerlimit;
    ebno=10^(k/10);
    snr(s,1)=k;
    for x=hc3/parts:hc3/parts:hc3-hc3/parts
        % zz=0.5*acos((2*x*(1/hc(aa)))-1);
      %  ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc3-x));                  % for normal distribution with narrow fov
      %  ff3=exp(((zz-(pi/6))^2)/(2*(pi/9)));
        ff=(2/pi)/(ff2);
        ggn=erfc(sqrt(ebno*x*0.5))*0.5*ff;
        %ggn=erfc(sqrt(ebno*x))*0.5*ff;
        pen(s,aaa)=pen(s,aaa)+ ggn*(hc3/parts);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

semilogy(snr,pen(:,1),'-+b');
hold on
semilogy(snr,pen(:,2),'-*g');
hold on
semilogy(snr,pen(:,3),'-or');
grid;
% title('BER.VS.SNR (Orientation of Receiver is Uniformly distributed)');
ylabel('BER');
xlabel('SNR [dB]');
legend('FOV=35 degree','FOV=60 degree','WIDE FOV');


