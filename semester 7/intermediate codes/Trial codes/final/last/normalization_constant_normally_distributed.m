clc
clear all;
p=0;
fov1=pi*(35/180);
fov2=pi*(60/180);
hc1= 6.7451e-13;
hc3= 6.7451e-13;
hc2= 6.7451e-13;
fov=[fov1 fov2 ];

product1=0.058161757735338 ;
cdf1=0.441287994957038;
product2= 1.849314141758325e-12;
cdf2= 0.187747628604906;
parts=1000000;

for x=hc3/parts:hc3/parts:hc3-(hc3/parts)
         zz=0.5*acos((2*x*(1/hc2))-1);
        ff1=sqrt(2*pi*(pi/9));
        ff2=sqrt((4*x)*(hc2-x));                  % for normal distribution with wide FOV
        ff3=exp(((zz-(pi/6))^2)/(2*(pi/9)));
        ff=(1)/(ff1*ff2*ff3);
        p=p+ ff*(hc3/parts);
end

    co3=1/p                                          % wide fov
    
    
    

    p=0;
    aa=1;
    ini=(hc2)*(cos(fov(aa))^2);
fin=hc2;
chang=(fin-ini)/parts;
  for x=ini+chang:chang:fin-chang
        f1=sqrt(2*pi*(pi/9));
        f2=sqrt((4*x)*(hc1-x));
        f3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        f=(1/(f1*f2*f3));
        p=p+ f*(chang);
  end  
% c01=(1-product1*(1-(cos(fov1))^2))/p                                % fov= 35 degree
c01=cdf1/p
  

  p=0;
    aa=2;
    ini=(hc2)*(cos(fov(aa))^2);
fin=hc2;
chang=(fin-ini)/parts;
  for x=ini+chang:chang:fin-chang
        f1=sqrt(2*pi*(pi/9));
        f2=sqrt((4*x)*(hc2-x));
        f3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        f=1/(f1*f2*f3);
        p=p+ f*(chang);
  end  
%  c02=(1-product2*(1-(cos(fov2))^2))/p                                  % fov= 60 degree
c02=cdf2/p
  
  
  
  