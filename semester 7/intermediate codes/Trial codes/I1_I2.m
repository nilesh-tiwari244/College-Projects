clc
clear all
syms l h t 
hc1=2*(10)^(-6);
hc=1;

a= 5*(h^2);
I1=0;
uu=0;
 for xx=0:1:50
    f= (t^(xx-1)) * exp(-t);
    III1=int(f,t,[0,a]);
    ee= (   (4^xx) / (factorial(xx)*gamma(xx+2)) ) ;
    jj = ( (  ee  * ( ((h^2)/0.1)^xx )   *     ( 1 / exp(40*(h^2) )  )  *  (  III1  )  ) );
    if (ee < 0.0001)
        uu=uu+1;
    end
    if (uu>2)
        break
    end
    I1= jj+I1;
end
abs(int(I1,h,0,1));

III2=int(f,t,[a,inf]);
ga2= (1/ factorial(1))* III2;
I2 = ga2;
fh=(2/(3.14*hc1))/ ((1- (h/hc1))^0.5);
p1=int((fh*I1),h,0,hc1)
