tic
%clc;
%close all;
%clear all;

check=0;


cubedim=10;
L=cubedim;
B=cubedim;                   
H=8;

x1=L/2;        % led position
y1=B/2; 

gamma=3; 
fov = [pi*(60/180) pi*(40/180)]; 
th=0;

position_sensi=0.25;
incre=position_sensi;
[x,y] = meshgrid(0:incre:L , 0:incre:B); 
mat=zeros(length(x));

%%      combined Probability distribution 


modex=7; % choose such that ((mode/cubedim)*10) is natural numbers
varx=13; %  natural no. between 2 and 13
modey=7;
vary=13;
 [distri_alpha_x, distri_beta_x]=func_5_beta_distri_parameter_march28trial_twoantenna_in_room(modex/cubedim,varx);
 toc 
 [distri_alpha_y, distri_beta_y]=func_5_beta_distri_parameter_march28trial_twoantenna_in_room(modey/cubedim,vary);
 toc

kk1=2;
indexx=2*position_sensi;
limit=length(x);
distri=zeros(limit,limit);

 while kk1<=limit
     kk2=2;
     indexy=2*position_sensi;
   while kk2<=limit
        distri(kk1,kk2)=(betacdf(indexx/L,distri_alpha_x,distri_beta_x)-betacdf((indexx-position_sensi)/L,distri_alpha_x,distri_beta_x)) * (betacdf(indexy/L,distri_alpha_y,distri_beta_y)-betacdf((indexy-position_sensi)/L,distri_alpha_y,distri_beta_y));
        indexy=indexy+position_sensi;
        kk2=kk2+1;
   end
    indexx=indexx+position_sensi;
    kk1=kk1+1;
 end
 distribution=(distri)/sum(sum(distri));
 distribution=distribution';
 
 %% x probability density
 
 kk1=2;
indexx=2*position_sensi;
limit=length(x);
distri=zeros(1,limit);

 while kk1<=limit
        distri(1,kk1)=(betacdf(indexx/L,distri_alpha_x,distri_beta_x)-betacdf((indexx-position_sensi)/L,distri_alpha_x,distri_beta_x)) ;
    indexx=indexx+position_sensi;
    kk1=kk1+1;
 end
 distribution_x=(distri)/sum(sum(distri));
 
 %% y probability density
distri_alpha_x=distri_alpha_y;
distri_beta_x=distri_beta_y;      
 kk1=2;
indexx=2*position_sensi;
limit=length(x);
distri=zeros(1,limit);

 while kk1<=limit
        distri(1,kk1)=(betacdf(indexx/L,distri_alpha_x,distri_beta_x)-betacdf((indexx-position_sensi)/L,distri_alpha_x,distri_beta_x)) ;
    indexx=indexx+position_sensi;
    kk1=kk1+1;
 end
 distribution_y=(distri)/sum(sum(distri));

%%
klo=0;
limii=10;
for ite=1:1:limii
datasi=1000;
angx1=99.3103   ;
angy1=93.1034;
angx2=74.4828  ;
angy2=   99.3103 ;
x=zeros(1,datasi);
y=zeros(1,datasi);
cheeee=0;
for l1=1:1:datasi
    x(1,l1)=func_6_generate_random(distribution_x,0,L);
    y(1,l1)=func_6_generate_random(distribution_y,0,B);
    cheeee=cheeee+1;
    if mod(cheeee,100)==0
        cheeee;
    end
end
[x,y] = meshgrid(x , y);  
a=x1-x;
b=y1-y;
c=H;
Sx1=cos(angx1*(pi/180));       
Sy1=cos(angy1*(pi/180));
Sx2=cos(angx2*(pi/180));       
Sy2=cos(angy2*(pi/180));
Sz1=sqrt(1-(Sx1^2 + Sy1^2));
Sz2=sqrt(1-(Sx2^2 + Sy2^2));
R=sqrt((a.^2)+(b.^2)+(c.^2));
psi=acos(H*(R.^(-1)));
theta1=acos(((Sx1*a)+(Sy1*b)+(Sz1*c)).*(R.^(-1)));
theta2=acos(((Sx2*a)+(Sy2*b)+(Sz2*c)).*(R.^(-1)));
h1=func_9_calculating_channel_coefficient(length(x),theta1,gamma,R,psi,fov(1));
h2=func_9_calculating_channel_coefficient(length(x),theta2,gamma,R,psi,fov(2));
hnet=((h1.^2)+(h2.^2)).^0.5;
exph2=sum(sum(hnet))/(datasi^2);
klo=klo+exph2;
ite
end
toc
klo/limii


%%                      Final Results
% h=    6.6139e-06;  % single antenna
% h= 7.39e-06;  % double antenna

%%     Single antenna

 figure
 func_8_simulation_ber_snr_plot(6.55e-06,90,'--ob');
 hold on
 func_7_analytical_ber_snr_plot(6.6139e-06,90,113,'-y');
 legend('Simulation','Analytical')
 %%    double antenna
 figure
  func_8_simulation_ber_snr_plot(7.35e-06,90,'--Xr');
 hold on
 func_7_analytical_ber_snr_plot(7.39e-06,90,113,'-g');
 legend('Simulation','Analytical');
%%     combined
 figure
 func_8_simulation_ber_snr_plot(6.55e-06,90,'--ob');
 hold on
 func_7_analytical_ber_snr_plot(6.6139e-06,90,114,'-y');
 hold on
 func_8_simulation_ber_snr_plot(7.35e-06,90,'--Xr');
 hold on
 func_7_analytical_ber_snr_plot(7.39e-06,90,114,'-g');
 legend('one photo-detector, Simulation','one photo-detector, Analytical','Two photo-detectors, Simulation','Two photo-detectors, Analytical');
%%


% [a1x,a1y] = meshgrid(0:incre:L , 0:incre:B);    
% surf(a1x,a1y,x.*y);
% xlabel('X coordinate');
% ylabel('Y coordinate');
% zlabel('Probability distribution'); 
    
    %%
    
% N=1000;
% width=60;
% lower_limit=60;
% incre=2.5;
% x=randi([0 1],1,N);
% mod_data=pskmod(x,2);
% BER=zeros(1,width/incre);
% new_data=zeros(1,N);
% 
% check_plot=0;
% 
% for SNR=lower_limit:incre:lower_limit+width
%     
% check_plot=check_plot+1;
% if mod(check_plot,10)==0
%         check_plot
%     end
% err1=0;
% en=1/(10^(SNR/10));
%  
% for tt=1:1:N
% 
% noi1=normrnd(0,sqrt(en));
% noi2=normrnd(0,sqrt(en));
% recent=mod_data(tt);
% new_data(1,:)= ((h1.*(noi1 + h1*(recent)))+ (h1.(noi2 + h2*(recent)))).* aas;
% 
% if (new_data(1,tt)>0) && (recent<0)
%     err1=err1+1;
% elseif (new_data(1,tt)<0) && (recent>0)
%        err1=err1+1;
%    else 
%        err1=err1+0;
% end
% end
% BER(1,check_plot)=(err1/N);
% end
% figure
% SNR=lower_limit:incre:lower_limit+width;
% ylabel('BER');
% xlabel('SNR [dB]');
% semilogy(SNR(:),BER(1,:),'-og');
% hold on
% grid on
    
    