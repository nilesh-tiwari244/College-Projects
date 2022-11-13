syms x;
p=0;
fov1=pi*(35/180);
fov2=pi*(60/180);
delta1=0.1318


     ff1=sqrt(2*pi*(pi/9));                    
     ff3=exp(((x-(pi/6))^2)/(2*(pi/9)));
     ff=1/(ff1*ff3);
     
     
     g1=double(1-(int(ff,x,cos(fov1),inf)));
     cdf1=1/g1     %    for Fov=35 degree
     
     
   
     g2=double(1-(int(ff,x,cos(fov2),inf)));
     cdf2=1/g2    %    for Fov=60 degree
     
     delta2=9.85e-12