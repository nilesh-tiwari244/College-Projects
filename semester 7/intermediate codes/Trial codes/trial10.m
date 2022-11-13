syms t
hc=2e-6;                                             
Eb=4;                                                            % defining variables and constants
M=2;                                                
snrdb=10;
No=Eb/(exp(snrdb/20));                      
fh=@(h)(3.1847e+05)/ sqrt(1- ((h/hc).^2));                      % due to motion of receiver
 a=@(h)(0.5/No)*(h.^2);
     f2= t * exp(-t);
    G2=@(h) int(f2,t,[a inf])*fh;      
    
    
      format long
   j4=integral(G2,0,hc);