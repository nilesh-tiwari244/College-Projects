clc
clear all;


snr_lower_limit=0;
snr_width=1;


pe=zeros(snr_width,3);
snrdb=zeros(snr_width,1);
P=zeros(snr_width,1);
hc=(2.0000e-06);
Eb=4;                                                            % defining variables and constants
M=2;
check1=0;

for k=1:1:snr_width
    snrdb(k,1)=snr_lower_limit+k;
    n0=10^(snrdb(k,1)/10);
    jj=0;
    jj2=0;
    acc_h=1e-8;
    
for h=(acc_h):(acc_h):(hc)-(acc_h) 
fh=(3.1847e+05)/ sqrt((1- (h*(5e+5))^2));
a= (0.5/n0)*(h^2);
I1=0;

for xx=1:1:100                                                                     % loop for summing up
ee = (Eb^xx) /  (  factorial(xx) * gamma(xx+M)  )    ;

%     f1= (t^(xx-1)) * exp(-t);                                                     % lower incomplete gamma function 
%     G1=int(f1,t,[0 a]); 

 G1=0;  
for t=a*(1e-3):a*(1e-3):a*(1-(a*(1e-3)))
    f1= (t^(xx-1)) * exp(-t);
    G1=G1+ f1*a*(1e-3);
end

    I1 =  I1 +   (   ee   *   ( ( (h^2)/n0 )^xx  )  *  (    1/exp( (4/n0)*(h^2) )   )   *  G1    );
end

j1=I1*fh;
jj=jj+ j1 * acc_h ;
 
end
pe(k,1)=jj ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for h=(acc_h):(acc_h):(hc)-(acc_h) 
fh=(3.1847e+05)/ sqrt((1- (h*(5e+5))^2));
a= (0.5/n0)*(h^2);
% f2= (t) * exp(-t=ll);                                                     % lower incomplete gamma function 
% G23=int(f2,ll,a,inf);
G23=exp(-a)*(a + 1);
I2=G23;
j2=I2*fh;
 jj2=jj2 + j2 * acc_h ;
end
G23
fh
pe(k,2)=jj2;
pe(k,3)=0.5*(pe(k,1)+pe(k,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

check1=check1+1
end
semilogy(snrdb(:,1),pe(:,3),'-+r');
hold on
semilogy(snrdb,pe(:,2),'-*g');
 hold on
 semilogy(snrdb,pe(:,3),'-ob');
 grid;
 legend('total','pe2','pe1');
pe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    
   
    
    
    
    