clc
clear all;
p=0;
fov1=pi*(35/180);
fov2=pi*(60/180);
hc1=(4e-12)/((sin(fov1))^2);
hc3=4e-12;
hc2=(4e-12)/((sin(fov2))^2);

product1=0.058161757735338 ;
cdf1=0.441287994957038;
product2= 1.849314141758325e-12;
cdf2= 0.187747628604906;

for x=hc3/2000:hc3/2000:hc3-(hc3/2000)
        f1=sqrt(2*pi*(pi/9));
        f2=sqrt((4*x)*(hc3-x));
        f3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        f=1/(f1*f2*f3);
        p=p+ f*(hc3/2000);
end
    co3=1/p                                          % wide fov
    
    
    

    p=0;
  for x=((cos(fov1))^2)*hc1:hc1/2000:hc1-(hc1/2000)
        f1=sqrt(2*pi*(pi/9));
        f2=sqrt((4*x)*(hc1-x));
        f3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        f=(1/(f1*f2*f3));
        p=p+ f*(hc1/2000);
  end  
% c01=(1-product1*(1-(cos(fov1))^2))/p                                % fov= 35 degree
c01=1/p
  

  p=0;
  for x=((cos(fov2))^2)*hc2:hc2/2000:hc2-(hc2/2000)
        f1=sqrt(2*pi*(pi/9));
        f2=sqrt((4*x)*(hc2-x));
        f3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        f=1/(f1*f2*f3);
        p=p+ f*(hc2/2000);
  end  
%  c02=(1-product2*(1-(cos(fov2))^2))/p                                  % fov= 60 degree
c02=1/p
  
  
  
  