clc
clear all
syms h t r
hc=2*(10)^(-6);
k1=((h^2)/0.1);
k2=(r/4)^(1/2);
k3=exp((r+4)*((h)^2)*10);
k4=20*((h)^2)*(4*r)^(1/2);
k5=(exp(k4*cos(t)))*(sin(t))^2;
k6=int(k5,t,[0,3.147]);
k7=k6*(k4/(2*3.147));
k8=k1*k2*(1/k3)*k7;
i1=int(k8,r,[0,0.5]);
s2=k1^2;
s3=1/exp((k1*r));
s4=s2*s3*r;
i2=int(s4,r,[0.5,inf]);
fh=(2/(3.14*hc))/ ((1- (h/hc))^0.5);
p1=int((fh*i1),h,0,inf)

