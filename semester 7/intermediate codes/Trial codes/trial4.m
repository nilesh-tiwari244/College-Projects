clc                                                 %  prerequisite               
syms h t
hc1=2.0000e-06;                                     % defining variables and constants
flag=0;
P1=0;                                               % prob. of error while sending 1
check=0;
%%
a= 5*(h^2);
fh=(3.1847e+05)/ sqrt((1- (h*(5e+5))^2))            % due to motion of receiver
%%
for xx=1:1:inf                                      % loop for summing up
    ee= ( (4^xx) / (factorial(xx)*gamma(xx+2)) ) ;
%%   
    if (ee < 0.01)
        flag=flag+1;
    end                                              % both these if loops are for breaking the infinite loop
    if (flag>2)
        break
    end
 %%
   check=check+1
   f= (t^(xx-1)) * exp(-t);
    G=int(f,t,[0 a]);                                % lower incomplete gamma function 
    I1 =((ee*(((h^2)/0.1)^xx)*(1/exp(40*(h^2) ))*(G)));
    j1=I1*fh;
    j2=int(j1,h,[0 hc1]);
    j3=vpa(j2);                                      % for getting numerical value from definite integral
    P1=P1+ j3;  
end
P1                                                   % final prob. of error when s=1 is send
I1


