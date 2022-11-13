tic

% clc;
% close all;
% clear all;


%%                                constant in channel coefficient expression
gamma=3;                      
Ar=1e-4;
g=1;

%%                                 dimension of room
cubedim=10;
L=10;
B=10;                   
H=8;

%%                                 FOV of antenna (Physical Limitation on antenna)
fov = pi*(60/180); 

%%                                 Position of Led at ceiling
x1=L/2;
y1=B/2;                

%%                                 Range of possible direction vector for antenna
angx=linspace(0,180,250);
angy=linspace(0,180,250);
ra1=cos(angx*(pi/180));       
ra2=cos(angy*(pi/180));

%%                                 VARIABLES for storing position of DIRECTION VECTOR LEADING TO maximum h
ttt=1;
cuma=0;
pos=[0 0];        
summ=zeros(1,length(ra1)^2);

%%                                  GRID OF Available positions to the receiver 
incre=0.1;    % decresing this smoothen the curves but also increase complexity hence time consumption
position_sensi=incre;
[x,y] = meshgrid(0:incre:L , 0:incre:B);    
h_max = zeros(length(x));
h = zeros(length(x)); 

%%                                 DEFINING PROBABILITY DISTRIBUTION FOR MOVEMENT OF THE RECEIVER  (BETA DISTRIBUTION)   
% distri_alpha_x=12;
% distri_beta_x=3;      %%% mode= (distri_alpha - 1 ) / (distri_alpha + distri_beta - 2)
% distri_alpha_y=2;
% distri_beta_y=8;

modex=7; % choose such that ((mode/cubedim)*L) is natural numbers
varx=13; %  natural no. between 2 and 13
modey=7;
vary=13;
 tic
 [distri_alpha_x, distri_beta_x]=func_5_beta_distri_parameter_march28trial_twoantenna_in_room(modex/cubedim,varx);
 toc 
 tic
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
 
%%                           FINDING POSITION VECTOR OF ANTENNA, LEADING TO MAXIMUM expected h 
check=0;
for qq=1:1:length(ra1)
     check=check+1;
        if mod(check,5)==0
          check
        end
    S1=ra1(qq);
    for jj=1:1:length(ra1)
        pp=jj;
        S2=ra2(pp);
        if (S2^2 + S1^2) > 1
           summ(1,ttt)=0;
        else
        S3=sqrt(1-(S1^2 + S2^2));
a=x1-x;
b=y1-y;
c=H;
R=sqrt((a.^2)+(b.^2)+(c.^2));
psi=acos(H*(R.^(-1)));
theta=acos(((S1*a)+(S2*b)+(S3*c)).*(R.^(-1)));
for i=1:1:length(x)
    for j=1:1:length(x)
if abs(theta(i,j)) < fov/2
    h(i,j)=((gamma+1)*(1/(R(i,j)^2))*(0.5*10e-4)*(1/pi)* cos(theta(i,j))*(cos(psi(i,j)))^(gamma));
else
    h(i,j)=0;
end
    end
end                                  

exph=h.*distribution;
summ(1,ttt)=sum(sum(exph));
        end
%%                          OPTIMISING S1,S2 AND S3 BY MAXIMISING EXPECTED h
        if cuma < summ(1,ttt)
        pos=[qq jj];
        cuma=summ(1,ttt);
        h_max=h;
         ttt=ttt+1;
        else   
         ttt=ttt+1;
        end 
    end
end
%%                             STORING S1,S2 AND S3 LEADING TO MAXIMUM EXPECTED h                             
S1max=ra1(pos(1));
S2max=ra2(pos(2));
S3max=sqrt(1-(S1max^2 + S2max^2));
anglex=acos(S1max)*(180/pi);
angley=acos(S2max)*(180/pi);
anglez=acos(S3max)*(180/pi);
S_max= [S1max S2max S3max];
angle_from_axes= [anglex angley anglez]

%%                             STORING MAXIMUM EXPECTED VALUE h 
expected_h=cuma
occ=zeros(length(x));
for l1=1:1:length(x)
    for l2=1:1:length(x)
        if(h_max(l1,l2)>0)
            occ(l1,l2)=1;
        end
    end
end
exp_coverage=sum(sum(occ.*distribution))
totalre=((cubedim/position_sensi)+1)^2;
sum_of_h=sum(sum(occ));
coverage=sum_of_h/totalre

%%                STORING COORDINATES OF RECEIVER, LEADING TO MAXIMUM h for best possible COMBINATION OF  S's
maxh=max(max(h_max));
[K1,K2]=find(h_max==maxh);
xcoordi=x(K1,K2);
ycoordi=y(K1,K2);
coords=[xcoordi ycoordi];

%%                                              PLOTTING GRAPHS
% figure
% mesh(x,y,distribution);
% xlabel('X coordinate');
% ylabel('Y coordinate');
% zlabel('Probability distribution');
figure
surf(x,y,h_max);
xlabel('X coordinate');
ylabel('Y coordinate');
zlabel('h values');

%%                                            PLOTTING THE CONE OF FOV

% figure
% func_2_cone_plotting_march28trial_twoantenna_in_room(fov,S1max, S2max, S3max,'y')
% func_8_simulation_ber_snr_plot(neth,60,'-or');
% toc
