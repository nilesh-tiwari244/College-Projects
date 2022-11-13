clc;
close all;

L=10;
B=10;                   % dimension of room
H=10;
gamma=3;

fov = pi*(60/180);    % FOV of antenna

x1=L/2;
y1=B/2;                % Position of Led at ceiling

incre=0.1;
[x,y] = meshgrid(0:incre:L , 0:incre:B);

anglex=90;
angley=90;

S1=cos(anglex*(pi/180));
S2=cos(angley*(pi/180));
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

 mesh(x,y,h);
 maxh=max(max(h));
 [K1,K2]=find(h==maxh);
 xcoordi=x(K1,K2);
 ycoordi=y(K1,K2);
 coords=[xcoordi ycoordi];
 sum_of_h=sum(sum(h))
 average_h=sum_of_h/(length(x)^2)
 angle_from_axes= [anglex angley acos(S3)*(180/pi)]

 xlabel('X coordinate');
 ylabel('Y coordinate');
 zlabel('h values');
       