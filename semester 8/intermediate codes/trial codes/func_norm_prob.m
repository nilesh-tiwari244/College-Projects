function [prob]= func_norm_prob(theta_mean,theta_var,low_bound,up_bound)

syms x;

ff1=sqrt(2*pi*(theta_var));                    
ff3=exp(((x-(theta_mean))^2)/(2*(theta_var)));
ff=1/(ff1*ff3);

prob=double((int(ff,x,low_bound,up_bound)));