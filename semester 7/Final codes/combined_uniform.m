clc
clear all;
close all;


dim=40;                                                     %%%% width of x axis
snrlowerlimit=105;                                           %%%% x axis is ber vs snr
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
gamma=3;                                 
d=2.5;
cdf=[cdf1 cdf2];
fov=[fov1 fov2 ];
const=[0.1334 0.3042 0.4057];


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
semilogy(snr,BER(:,2),'-b');
hold on
semilogy(snr,BER(:,3),'-r');
grid on
hold on

ylabel('BER');
xlabel('SNR [dB]');
% legend('FOV=35 degree','FOV=60 degree','WIDE FOV');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;
clear all;
%close all;

N=10000;
width=40;
lower_limit=105;
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
semilogy(SNR(:),BER(1,:),'--og');
hold on
semilogy(SNR(:),BER(2,:),'--+b');
hold on;
semilogy(SNR(:),BER(3,:),'--*r');
legend('FOV=35 degree(Analytical)','FOV=60 degree(Analytical)','WIDE FOV(Analytical)','FOV=35 degree(Simulation)','FOV=60 degree(Simulation)','WIDE FOV(Simulation)');

legend('Location','southwest')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

