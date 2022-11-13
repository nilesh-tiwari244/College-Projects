run=0;
p=0;
l=4e-12;
for x=l/2000:l/2000:l/2
        f1=sqrt(2*pi*(pi/9));
        f2=sqrt((4*x)*(l-x));
        f3=exp(((x-(pi/6))^2)/(2*(pi/9)));
        f=1/(f1*f2*f3);
        run=run+1;
        p=p+ f*(l/2000);
end
    k=1/p