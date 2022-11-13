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
angx1=93.9759   ;
angy1=93.9759;
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
Sz1=sqrt(1-(Sx1^2 + Sy1^2));
R=sqrt((a.^2)+(b.^2)+(c.^2));
psi=acos(H*(R.^(-1)));
theta1=acos(((Sx1*a)+(Sy1*b)+(Sz1*c)).*(R.^(-1)));
h1=func_partial_effici_march28trial_twoantenna_in_room(length(x),theta1,gamma,R,psi,fov(1));
hnet=((h1));
exph1=sum(sum(hnet))/(datasi^2);
klo=klo+exph1;
ite
end
toc
klo/limii

%  figure
%  func_plot_ber_vs_snr(exph,60,'--*r');
%  hold on
%  semilogy(snrana(:,1),BERana(:,1),'-g');
%  legend('Simulation','Analytical');