function []=func_effici_march28trial_twoantenna_in_room(qq1,jj1,qq2,jj2,Sx1, Sy1, Sx2, Sy2)

global distribution Ar g threshold L B H angle_sensi position_sensi x1 y1 angx1 angy1 angx2 angy2 rax1 ray1 rax2 ray2;
global ttt cuma pos summ x y h_max incre gamma fov ;

 if (Sy1^2 + Sx1^2) > 1
           summ(1,ttt)=0;
        elseif (Sy2^2 + Sx2^2) > 1
           summ(1,ttt)=0;
        else
        Sz1=sqrt(1-(Sx1^2 + Sy1^2));
        Sz2=sqrt(1-(Sx2^2 + Sy2^2));
h1 = zeros(length(x)); 
h2 = zeros(length(x));
a=x1-x;
b=y1-y;
c=H;
R=sqrt((a.^2)+(b.^2)+(c.^2));
psi=acos(H*(R.^(-1)));
theta1=acos(((Sx1*a)+(Sy1*b)+(Sz1*c)).*(R.^(-1)));
theta2=acos(((Sx2*a)+(Sy2*b)+(Sz2*c)).*(R.^(-1)));


h1=func_partial_effici_march28trial_twoantenna_in_room(length(x),theta1,gamma,R,psi,fov);
h1=func_partial_effici_march28trial_twoantenna_in_room(length(x),theta2,gamma,R,psi,fov);

% for i=1:1:length(x)
%     for j=1:1:length(x)
% if abs(theta1(i,j)) < fov/2
%     h1(i,j)=((gamma+1)*(1/(R(i,j)^2))*(0.5*10e-4)*(1/pi)* cos(theta1(i,j))*(cos(psi(i,j)))^(gamma));
% else
%     h1(i,j)=0;
% end
%     end
% end
% 
% for i=1:1:length(x)
%     for j=1:1:length(x)
% if abs(theta1(i,j)) < fov/2
%     h2(i,j)=((gamma+1)*(1/(R(i,j)^2))*(0.5*10e-4)*(1/pi)* cos(theta2(i,j))*(cos(psi(i,j)))^(gamma));
% else
%     h2(i,j)=0;
% end
%     end
% end
                        %OPTIMISING S1,S2 AND S3 BY MAXIMISING EXPECTED h
[summ(1,ttt),h1]=abov_thre(h1,h2,length(x),length(x),threshold,distribution);
        end
        if cuma < summ(1,ttt)
        pos=[qq1 jj1 qq2 jj2];
        cuma=summ(1,ttt);
        h_max=h1;
         ttt=ttt+1;
        else   
         ttt=ttt+1;
        end 