clc                                                             %  prerequisite               
syms  t 
hc=2e-6;                                             
Eb=4;                                                            % defining variables and constants
M=2;                                                
snrdb=10;
No=Eb/(exp(snrdb/20));
tt=0;
ll=0;
for h=(0.0001e-6):1e-8:(1e-6)                      
      fh=(3.1847e+05)/ sqrt((1- (h/(hc))^2));                      % due to motion of receiver
   
  
    ll=ll+1
    a=double((0.5/No)*(h^2));
    
    
    f2= t * exp(-t);
    G2=double(int(f2,t,[a inf]));
    I2= double(G2/factorial(M-1));
    
    j2=double(I2*fh);
  
    tt=tt+(j2*(1e-8));
end
tt
    
   