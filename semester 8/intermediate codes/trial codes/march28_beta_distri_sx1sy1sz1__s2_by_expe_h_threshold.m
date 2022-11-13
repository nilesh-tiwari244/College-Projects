tic

clc;
close all;
clear all;


%%                                constant in channel coefficient expression
gamma=3;                      
Ar=1e-4;
g=1;

%%                                 dimension of room
L=10;
B=10;                   
H=10;

%%                                 FOV of antenna (Physical Limitation on antenna)
fov = pi*(50/180); 

%%                                 Position of Led at ceiling
x1=L/2;
y1=B/2;                

%%                                 Range of possible direction vector for antenna
angx=linspace(0,180,181);
angy=linspace(0,180,181);
ra1=cos(angx*(pi/180));       
ra2=cos(angy*(pi/180));

%%                                 VARIABLES for storing position of DIRECTION VECTOR LEADING TO maximum h
ttt=1;
cuma=0;
pos=[0 0];        
summ=zeros(1,length(ra1)^2);

%%                                  GRID OF Available positions to the receiver 
incre=0.1;    % decresing this smoothen the curves but also increase complexity hence time consumption
[x,y] = meshgrid(0:incre:L , 0:incre:B);    
h_max = zeros(length(x));

%%                                 DEFINING PROBABILITY DISTRIBUTION FOR MOVEMENT OF THE RECEIVER  (BETA DISTRIBUTION)   
distri_alpha_x=12;
distri_beta_x=3;      %%% mode= (distri_alpha - 1 ) / (distri_alpha + distri_beta - 2)
distri_alpha_y=2;
distri_beta_y=8;

kk1=2;
indexx=2;
indexy=2;
limit=length(x);
distri=zeros(length(x),length(x));

 while kk1<=limit
     kk2=2;
     indexy=2;
   while kk2<=limit
        distri(indexx,indexy)=(betacdf((1/L)*x(1,indexx),distri_alpha_x,distri_beta_x)-betacdf((1/L)*x(1,indexx-1),distri_alpha_x,distri_beta_x)) * (betacdf((1/B)*y(indexy,1),distri_alpha_y,distri_beta_y)-betacdf((1/B)*y(indexy-1,1),distri_alpha_y,distri_beta_y));
        indexy=indexy+1;
        kk2=kk2+1;
   end
    indexx=indexx+1;
    kk1=kk1+1;
 end
 sum(sum(distri));
 distribution=(distri)/sum(sum(distri));
 
 
%%                           FINDING POSITION VECTOR OF ANTENNA, LEADING TO MAXIMUM expected h 
for qq=1:1:length(ra1)
    S1=ra1(qq);
    for jj=1:1:length(ra1)
        pp=jj;
        S2=ra2(pp);
        if (S2^2 + S1^2) > 1
           summ(1,ttt)=0;
        else
        S3=sqrt(1-(S1^2 + S2^2));
h = zeros(length(x)); 
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
%%                          OPTIMISING S1,S2 AND S3 BY MAXIMISING EXPECTED h
exph=h.*distribution;
 summ(1,ttt)=sum(sum(exph));
        end
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
expected_h=cuma;
sum_of_h=sum(sum(h_max));
average_h=sum_of_h/(length(x)^2);

%%                STORING COORDINATES OF RECEIVER, LEADING TO MAXIMUM h for best possible COMBINATION OF  S's
maxh=max(max(h_max));
[K1,K2]=find(h_max==maxh);
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

%%                                            PLOTTING THE CONE OF FOV

rad=tan((fov/2));
figure
p1=[0 0 0];
p2=1.5*[S1max S2max S3max];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'k')
hold on;
axis_length=1.5;
p2=[0 0 axis_length];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'b')
hold on;
p2=[0 axis_length 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'g')
hold on;
p2=[axis_length 0 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'r')
hold on;
%% 
p2=[0 0 -1*axis_length];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'b')
hold on;
p2=[0 -1*axis_length 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'g')
hold on;
p2=[-1*axis_length 0 0];
v=[p2;p1];
plot3(v(:,1),v(:,2),v(:,3),'r')
hold on;
pre=0.01;
for xx=-2:pre:2
    for yy=-2:pre:2
        for zz=-2:pre:2
            pp=S1max*(xx-S1max)+S2max*(yy-S2max)+S3max*(zz-S3max);
            if abs(pp)<0.01
                dis=sqrt((xx-S1max)^2+(yy-S2max)^2+(zz-S3max)^2);
                if abs(dis-rad)<0.01
                    p2=[xx yy zz];
                    v=[p2;p1];
                    plot3(v(:,1),v(:,2),v(:,3),'y')
                    hold on;
                end
            end
        end
    end
end
legend('Direction vector','Z axis','Y axis','X axis');
xlabel('X coordinate');
ylabel('Y coordinate');
zlabel('Z coordinate');

toc
