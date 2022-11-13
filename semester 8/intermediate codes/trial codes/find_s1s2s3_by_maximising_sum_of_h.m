clc;
close all;

L=10;
B=10;                   % dimension of room
H=10;
gamma=3;

fov = pi*(60/180);    % FOV of antenna

x1=L/2;
y1=B/2;                % Position of Led at ceiling


angx=linspace(0,180,181);
angy=linspace(0,180,181);
ra1=cos(angx*(pi/180));
ra2=cos(angy*(pi/180));
summ=zeros(1,length(ra1)^2);

ttt=1;
cuma=0;
pos=[0 0];

incre=1;
[x,y] = meshgrid(0:incre:L , 0:incre:B);
h_max = zeros(length(x));

for qq=1:1:length(ra1)
    S1=ra1(qq);
    for jj=1:1:length(ra1)
        pp=jj;
        S2=ra2(pp);
        if (S2^2 + S1^2) > 1
           summ(ttt)=0;
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
 summ(1,ttt)=sum(sum(h));
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
S1max=ra1(pos(1));
S2max=ra2(pos(2));
S3max=sqrt(1-(S1max^2 + S2max^2));
anglex=acos(S1max)*(180/pi);
angley=acos(S2max)*(180/pi);
anglez=acos(S3max)*(180/pi);
angle_from_axes= [anglex angley anglez]
max_sum=cuma
 maxh=max(max(h_max));
 [K1,K2]=find(h_max==maxh);
 xcoordi=x(K1,K2);
 ycoordi=y(K1,K2);
 coords=[xcoordi ycoordi]
mesh(x,y,h_max);
xlabel('X coordinate');
ylabel('Y coordinate');
zlabel('h values');
