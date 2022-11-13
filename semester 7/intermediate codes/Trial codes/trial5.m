

clc                                                             %  prerequisite               
syms h t


hc=(2.0000e-06);                                               
flag=0;                                              
check=0;
Eb=4;                                                            % defining variables and constants
M=2;
snrdb=zeros(10,1);
P=zeros(10,1);                                                 

for h=(0.0001e-6):1e-8:(1e-6) 
fh=(3.1847e+05)/ sqrt((1- (h*(5e+5))^2));                        % due to motion of receiver
for loop=1:1:10
    
    snrdb(loop)=loop;
    Pe1=0;
    Pe2=0;
    No=Eb/(exp(loop/20));
    a= (0.5/No)*(h^2);
    f2= t * exp(-t);
    G2=int(f2,t,[a inf]);
    I2= G2/factorial(M-1);
    j2=I2*fh;
    jj2=int(jj1,h,[0 hc]);
    jjj2=vpa(jj2);                                             
    Pe2=Pe2+jjj2;
    
for xx=1:1:inf                                                                     % loop for summing up
    
    ee =   (4^xx) /   (  factorial(xx) * gamma(xx+M)  )    ;
  
    if (ee < 0.01)
        flag=flag+1;
    end                                                                            % both these if loops are for breaking the infinite loop
    if (flag>2)
        break
    end
    
    f1= (t^(xx-1)) * exp(-t);                                                     % lower incomplete gamma function 
    G1=int(f1,t,[0 a]);  
    
    I1 =     (   ee   *   ( ( (h^2)/No )^xx  )  *  (    1/exp( (4/No)*(h^2) )   )   *  G1    );
    j1=I1*fh;
    jj1=int(j1,h,[0 hc]);
    jjj1=vpa(jj1);                                                                 % for getting numerical value from definite integral
    Pe1=Pe1+jjj1; 
       
end
 
    P(loop,1)=0.5*abs(Pe1+Pe2); 

end
semilogy(snrdb,P,'-+g')

