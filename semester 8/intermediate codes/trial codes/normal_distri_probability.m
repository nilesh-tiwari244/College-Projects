function []= func_norm_prob(theta_mean,theta_var,low_bound,up_bound)

syms x;

theta_mean=pi*(25/180);
theta_var=0.01;

     ff1=sqrt(2*pi*(theta_var));                    
     ff3=exp(((x-(theta_mean))^2)/(2*(theta_var)));
     ff=1/(ff1*ff3);

     fov1=pi*(40/180);
     cdf1=double((int(ff,x,(fov1),inf)))