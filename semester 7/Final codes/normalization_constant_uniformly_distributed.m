clc
clear all;
close all;

fov1=pi*(35/180);
fov2=pi*(60/180);

fov=[fov1 fov2 ];

%%% uniform between 0 and 90 degrees                       
%%%% range of x=0 to 1

cdf1=0.6111;                                                   %%% 3----wide FOV
cdf2= 0.3333;

cdf=[cdf1 cdf2];
parts=1000000;

p=0;
for x=1/parts:1/parts:1-(1/parts)
        ff1=sqrt((4*x)*(1-x));                  % for uniform distribution with wide FOV
        ff2=2/pi;
        ff=(1)/(ff1*ff2);
        p=p+ ff*(1/parts);
end
    co3=1/p                                          % wide fov
   
    
    
    p=0;
    aa=1;
    ini=(cos(fov(aa))^2);
fin=1;
chang=(fin-ini)/parts;
  for x=ini+chang:chang:fin-chang
         ff1=sqrt((4*x)*(1-x));                  % for uniform distribution with narrow FOV
        ff2=2/pi;
        ff=(1)/(ff1*ff2);
        p=p+ ff*(1/parts);
  end  
                              % fov= 35 degree
c01=(1-cdf1)/(p)
  


  p=0;
    aa=2;
    ini=(cos(fov(aa))^2);
fin=1;
chang=(fin-ini)/parts;
  for x=ini+chang:chang:fin-chang
         ff1=sqrt((4*x)*(1-x));                  % for uniform distribution with narrow FOV
        ff2=2/pi;
        ff=(1)/(ff1*ff2);
        p=p+ ff*(1/parts);
  end  
                              % fov= 60 degree
c02=(1-cdf2)/(p)
  
  
  
  