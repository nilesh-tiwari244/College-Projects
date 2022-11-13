tic

clc;
close all;
clear all;

% global distribution Ar g threshold L B H angle_sensi position_sensi x1 y1 angx1 angy1 angx2 angy2 rax1 ray1 rax2 ray2;
% global ttt cuma pos summ x y h_max incre gamma fov ;

%%                                constant in channel coefficient expression
gamma=3;                      
Ar=1e-4;
g=1;
threshold=2e-6;

%%                                 dimension of room
cubedim=3;
L=cubedim;
B=cubedim;                   
H=2.4;
x1=L/2;        % led position
y1=B/2;      

%%                                 FOV of antenna (Physical Limitation on antenna)
fov = [pi*(40/180) pi*(40/180)]; 

%%                                 SELECT SENSITIVITY LEVEL
angle_sensi=10;
position_sensi=0.1;

% no_loops=((angle_sensi^4)+ (angle_sensi^2))/2;
% expectedtime=(0.5)*(angle_sensi/10)^4;


%%                                 Range of possible direction vector for antenna
angx1=linspace(0,180,angle_sensi);
angy1=linspace(0,180,angle_sensi);
rax1=cos(angx1*(pi/180));       
ray1=cos(angy1*(pi/180));
angx2=angx1;
angy2=angy1;
rax2=rax1;       
ray2=ray1;
%%                                 VARIABLES for storing position of DIRECTION VECTOR LEADING TO maximum h
ttt=1;
cuma=0;
pos=[1 1 1 1];        
summ=zeros(1,length(rax1)^4);
neth=0;
trialh=0;

%%                                  GRID OF Available positions to the receiver 
incre=position_sensi;    
[x,y] = meshgrid(0:incre:L , 0:incre:B);    
h_max = zeros(length(x));

%%                                 DEFINING PROBABILITY DISTRIBUTION FOR MOVEMENT OF THE RECEIVER  (BETA DISTRIBUTION)  

% distri_alpha_x=12;
% distri_beta_x=3;      %%% mode= (distri_alpha - 1 ) / (distri_alpha + distri_beta - 2)
% distri_alpha_y=2;
% distri_beta_y=8;

modex=7; % choose such that ((mode/cubedim)*10) is natural numbers
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
tic
for qq1=1:1:length(rax1)
    Sx1=rax1(qq1);
    
        check=check+1;
        if(cuma==1)
            break
        end
    if mod(check,2)==0
        check
    end
    
    for jj1=1:1:length(rax1)
        Sy1=ray1(jj1);
         if(cuma==1)
            break
         end
  
        
           for qq2=1:1:length(rax2)
        %for qq2=qq1:1:length(rax2)
        Sx2=rax2(qq2);
         if(cuma==1)
            break
        end
    
            if (Sy1^2 + Sx1^2) > 1
               for nax=1:1:length(rax2)
                  for pax=1:1:length(rax2)
                      summ(1,ttt)=0;
                       ttt=ttt+ 1;
                  end
               end
                break
            end
            
        if (qq1==qq2)
            st=jj1;
        else
            st=1;
        end
        
        
                 for jj2=1:1:length(rax2)
               % for jj2=st:1:length(rax2)
                      Sy2=ray1(jj2);
                       if(cuma==1)
                      break
                         end
                                
        if (Sy1^2 + Sx1^2) > 1
           summ(1,ttt)=0;
        elseif (Sy2^2 + Sx2^2) > 1
           summ(1,ttt)=0;
        else
        Sz1=sqrt(1-(Sx1^2 + Sy1^2));
        Sz2=sqrt(1-(Sx2^2 + Sy2^2));
h1 = zeros(length(x)); 
h2 = zeros(length(x));
a=x1-x;
b=y1-y;
c=H;
R=sqrt((a.^2)+(b.^2)+(c.^2));
psi=acos(H*(R.^(-1)));
theta1=acos(((Sx1*a)+(Sy1*b)+(Sz1*c)).*(R.^(-1)));
theta2=acos(((Sx2*a)+(Sy2*b)+(Sz2*c)).*(R.^(-1)));


h1=func_partial_effici_march28trial_twoantenna_in_room(length(x),theta1,gamma,R,psi,fov(1));
h2=func_partial_effici_march28trial_twoantenna_in_room(length(x),theta2,gamma,R,psi,fov(2));

                     %OPTIMISING S1,S2 AND S3 BY MAXIMISING EXPECTED h
[summ(1,ttt),h1,trialh]=abov_thre(h1,h2,length(x),length(x),threshold,distribution);
        end
        if cuma < summ(1,ttt)
        pos=[qq1 jj1 qq2 jj2];
        cuma=summ(1,ttt);
         ttt=ttt+1;
         h_max=h1;
         neth=trialh;
        else   
         ttt=ttt+1;
        end 
       
       
       % func_1_max_s1_s2_march28trial_twoantenna_in_room(qq1,jj1,qq2,jj2,Sx1, Sy1, Sx2, Sy2);
            end
        end
    end
end
toc
neth
%%                             STORING S1,S2 AND S3 LEADING TO MAXIMUM EXPECTED h                             
Sx1max=rax1(pos(1));
Sy1max=ray1(pos(2));
Sz1max=sqrt(1-(Sx1max^2 + Sy1max^2));
anglex1=acos(Sx1max)*(180/pi);
angley1=acos(Sy1max)*(180/pi);
anglez1=acos(Sz1max)*(180/pi);
S_max1= [Sx1max Sy1max Sz1max];
angle_from_axes1= [anglex1 angley1 anglez1];

Sx2max=rax2(pos(3));
Sy2max=ray2(pos(4));
Sz2max=sqrt(1-(Sx2max^2 + Sy2max^2));
anglex2=acos(Sx2max)*(180/pi);
angley2=acos(Sy2max)*(180/pi);
anglez2=acos(Sz2max)*(180/pi);
S_max2= [Sx2max Sy2max Sz2max];
angle_from_axes2= [anglex2 angley2 anglez2];
%%                             STORING MAXIMUM EXPECTED VALUE h 
expected_h=cuma
totalre=((cubedim/position_sensi)+1)^2;
sum_of_h=sum(sum(h_max));
 coverage=sum_of_h/totalre

%%                STORING COORDINATES OF RECEIVER, LEADING TO MAXIMUM h for best possible COMBINATION OF  S's
 maxh=max(max(distribution)); 
[K1,K2]=find(distribution==maxh);
xcoordi=x(K1,K2);
ycoordi=y(K1,K2);
coords=[xcoordi ycoordi];

%%                                              PLOTTING GRAPHS
mesh(x,y,distribution);
xlabel('X coordinate');
ylabel('Y coordinate');
zlabel('Probability distribution');
figure
mesh(x,y,h_max);
xlabel('X coordinate');
ylabel('Y coordinate');
zlabel('h values');
toc

%%                                            PLOTTING THE CONE OF FOV
% tic
% figure
% func_2_cone_plotting_march28trial_twoantenna_in_room(fov(1),Sx1max, Sy1max, Sz1max,'y')
% hold on
% func_2_cone_plotting_march28trial_twoantenna_in_room(fov(2),Sx2max, Sy2max, Sz2max,'c')
% toc