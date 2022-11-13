clc
clear all;
dim=45; 
snrlowerlimit=90;                        
                                         %%%% 1 is 35    
                                         %%%% 2 is 60
                                         %%%% 3 is wide 


fov1=pi*(35/180);
fov2=pi*(60/180);
fov3=pi/2;

pen=zeros(dim,3);
snr=zeros(dim,1);                                         



product1=0.058161757735338 ;
cdf1=0.441287994957038;
product2= 1.849314141758325e-12;
cdf2= 0.187747628604906;



fov=[fov1 fov2 fov3];
cdf=[cdf1 cdf2];
product=[product1 product2 0];
const=[3.68 2.12 1.42];




hc=[(4e-12)/((sin(fov1))^2) (4e-12)/((sin(fov1))^2) 4e-12];
hc3=4e-12;



for aa=1:1:2
    
for s=1:1:dim
    k=s+snrlowerlimit;
    en=10^(k/10);
    snr(s,1)=k;
    for x=hc(aa)*((cos(fov(aa))^2)+(1/2000)):hc(aa)/2000:hc(aa)-hc(aa)/2000
        
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc(aa)-x));                  % for normal distribution
        ff3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        ff=product(aa)/hc(aa)+(const(aa))/(ff1*ff2*ff3);
        ggn=erfc(sqrt(en*x*0.5))*0.5*ff;
        pen(s,aa)=pen(s,aa)+ ggn*(hc(aa)/2000);
        
    end
%      for x=hc(aa)/2000:hc(aa)/2000:hc(aa)*((cos(fov(aa))^2))  
%           
%         ffff=product(aa)/hc(aa);
%         ggggn=erfc(sqrt(en*x/2))*0.5*ffff;
%         pen(s,aa)=pen(s,aa)+ ggggn*(hc(aa)/2000);
%         
%     end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aaa=3;
for s=1:1:dim
    k=s+snrlowerlimit;
    en=10^(k/10);
    snr(s,1)=k;
    for x=hc3/2000:hc3/2000:hc3
        
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc3-x));                  % for normal distribution
        ff3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        ff=product(aaa)/hc3+(const(aaa))/(ff1*ff2*ff3);
        ggn=erfc(sqrt(en*x*0.5))*0.5*ff;
        pen(s,aaa)=pen(s,aaa)+ ggn*(hc3/2000);
        
    end
    
    
end

semilogy(snr,pen(:,1),'-or');
hold on
semilogy(snr,pen(:,2),'-*g');
hold on
semilogy(snr,pen(:,3),'-+b');
grid;
title('BER.VS.SNR (motion of receiver is normally distributed)');
ylabel('BER');
xlabel('SNR [dB]');
legend('FOV=35 degree','FOV=60 degree','wide FOV');


