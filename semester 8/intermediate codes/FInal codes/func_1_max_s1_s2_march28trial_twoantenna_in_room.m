function []=func_1_max_s1_s2_march28trial_twoantenna_in_room(qq1,jj1,qq2,jj2,Sx1, Sy1, Sx2, Sy2)

global distribution threshold H x1 y1 ;
global ttt cuma pos summ x y h_max gamma fov ;

 if (Sy2^2 + Sx2^2) > 1
  summ(1,ttt)=0;
  ttt=ttt+1;
 else
 Sz1=sqrt(1-(Sx1^2 + Sy1^2));
 Sz2=sqrt(1-(Sx2^2 + Sy2^2));
h1 = zeros(length(x)); 
h2 = zeros(length(x));
hresul=zeros(length(x));
a=x1-x;
b=y1-y;
c=H;
R=sqrt((a.^2)+(b.^2)+(c.^2));
psi=acos(H*(R.^(-1)));
theta1=acos(((Sx1*a)+(Sy1*b)+(Sz1*c)).*(R.^(-1)));
theta2=acos(((Sx2*a)+(Sy2*b)+(Sz2*c)).*(R.^(-1)));
h1=func_4_channel_coeff_march28trial_twoantenna_in_room(length(x),theta1,gamma,R,psi,fov);
h2=func_4_channel_coeff_march28trial_twoantenna_in_room(length(x),theta2,gamma,R,psi,fov);
[summ(1,ttt),hresul]=func_3_above_thre_march28trial_twoantenna_in_room(h1,h2,length(x),length(x),threshold,distribution);

        if cuma < summ(1,ttt)
        pos=[qq1 jj1 qq2 jj2];
        cuma=summ(1,ttt);
        h_max=hresul;
         ttt=ttt+1;
        else   
         ttt=ttt+1;
        end 
end
