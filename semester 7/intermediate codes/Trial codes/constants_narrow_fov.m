clc
clear all;
syms x;

     ff1=sqrt(2*pi*(pi/9));                    
     ff3=exp(((x-(pi/6))^2)/(2*(pi/9)));
     ff=1/(ff1*ff3);
     
                    %    for Fov=35 degree
     fov1=pi*(35/180)
     delta1=0.1318
     cdf1=double((int(ff,x,(fov1),inf)))    
     product1=(1-cdf1)*delta1
     
                    %    for Fov=60 degree
     fov2=pi*(60/180)               
     delta2=9.85e-12              
     cdf2=double((int(ff,x,(fov2),inf)))
     product2=(1-cdf2)*delta2
     
     
     
     