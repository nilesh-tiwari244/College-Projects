
clc;
clear all;
%close all;
parts=1000;
l=3;
limit=l^(-6);

snrlowerlimit=60;
dim=80;
incre=5;

hc2= 4e-12;
cd=0.1666;

BER=zeros(dim/incre,1);
snr=zeros(dim/incre,1);
check0=0;
for s=incre:incre:dim
    check0=check0+1;
    k=s+snrlowerlimit;
    en=10^(k/10);
    snr(check0,1)=k;
    penn=0;
    
    check3=0;
    check2=0;
    % *cos(pi*(35/180))   hc2*limit*(1-1/parts)
    for x=(hc2/parts)*(limit/parts):(hc2/parts)*(limit/parts):(hc2)*(limit/parts)*(1-(1/parts))
        g2=0;
        check2=check2+1;
        for y=limit/parts:limit/parts:limit-(limit/parts)
        check3=check3+1;
        check1=check2+(check3)*i;
        z=x/y;
        aunt=(2*z*(1/hc2));
        zz=0.5*acos(aunt-1);
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*z)*(hc2-z));                  % for normal distribution with wide FOV
        ff3=exp(((zz-(pi/6))^2)/(2*(pi/9)));
        ff=3.7/(ff1*ff2*ff3);
        part1=sqrt((y^(-0.33333))-(l^2));
        fd=(y^(-1.3333)) * (1/part1) *  (part1/4) * (1/exp((part1^2)/8)) ;
        g1=(1/y)*(cd*fd)* ff;
        g2=g2+g1*(limit/parts);
        end
        ggn=erfc(sqrt(en*x*0.5))*0.5*g2;
        penn=penn + ggn*(hc2/parts)*(limit/parts);
    end
        x=0; 
        ffff=0.44;
        ggggn=erfc(sqrt(en*x*0.5))*0.5*ffff;
        penn=penn+ ggggn; 
        BER(check0,1)=penn;
end
semilogy(snr(:,1),BER(1,:),'-+b');
BER