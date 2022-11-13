tic
clc;
close all;
clear all;

check=0;
cubedim=10;
L=cubedim;
B=cubedim;                   
H=8;
x1=L/2;        % led position
y1=B/2; 

gamma=3; 
fov = [pi*(60/180) pi*(30/180)]; 
th=0;

% syms x y ;
syms angx1 angy1 angx2 angy2; 
position_sensi=1;
incre=position_sensi;
[x,y] = meshgrid(0:incre:L , 0:incre:B); 
mat=zeros(length(x));

%%      combined Probability distribution 


modex=7; % choose such that ((mode/cubedim)*10) is natural numbers
varx=13; %  natural no. between 2 and 13
modey=7;
vary=13;
 tic
 [distri_alpha_x, distri_beta_x]=func_5_beta_distri_parameter_march28trial_twoantenna_in_room(modex/cubedim,varx);
 toc 
 tic
 [distri_alpha_y, distri_beta_y]=func_5_beta_distri_parameter_march28trial_twoantenna_in_room(modey/cubedim,vary);
 toc

kk1=2;
indexx=2*position_sensi;
limit=length(x);
distri=zeros(limit,limit);

 while kk1<=limit
     kk2=2;
     indexy=2*position_sensi;
   while kk2<=limit
        distri(kk1,kk2)=(betacdf(indexx/L,distri_alpha_x,distri_beta_x)-betacdf((indexx-position_sensi)/L,distri_alpha_x,distri_beta_x)) * (betacdf(indexy/L,distri_alpha_y,distri_beta_y)-betacdf((indexy-position_sensi)/L,distri_alpha_y,distri_beta_y));
        indexy=indexy+position_sensi;
        kk2=kk2+1;
   end
    indexx=indexx+position_sensi;
    kk1=kk1+1;
 end
 distribution=(distri)/sum(sum(distri));
 distribution=distribution';
 
 %% x probability density
 
 kk1=2;
indexx=2*position_sensi;
limit=length(x);
distri=zeros(1,limit);

 while kk1<=limit
        distri(1,kk1)=(betacdf(indexx/L,distri_alpha_x,distri_beta_x)-betacdf((indexx-position_sensi)/L,distri_alpha_x,distri_beta_x)) ;
    indexx=indexx+position_sensi;
    kk1=kk1+1;
 end
 distribution_x=(distri)/sum(sum(distri));
 
 %% y probability density
distri_alpha_x=distri_alpha_y;
distri_beta_x=distri_beta_y;      
 kk1=2;
indexx=2*position_sensi;
limit=length(x);
distri=zeros(1,limit);

 while kk1<=limit
        distri(1,kk1)=(betacdf(indexx/L,distri_alpha_x,distri_beta_x)-betacdf((indexx-position_sensi)/L,distri_alpha_x,distri_beta_x)) ;
    indexx=indexx+position_sensi;
    kk1=kk1+1;
 end
 distribution_y=(distri)/sum(sum(distri));
%  che=[0.1 0.2 0.3 0.4];
%  kkk=func_6_generate_random(che)

 %%

% distri_alpha_x = 3.0769 ;
% distri_beta_x= 1.8846;
% distri_alpha_y = 3.0769 ;
% distri_beta_y= 1.8846;
% 
% px=x/L;
% py=y/B;
% pdfx=vpa(((px)^(distri_alpha_x-1))*((1-px)^(distri_beta_x-1))*gamma(distri_alpha_x+distri_beta_x)*(1/gamma(distri_beta_x))*(1/gamma(distri_alpha_x)),4);
% pdfy=vpa(((py)^(distri_alpha_y-1))*((1-py)^(distri_beta_y-1))*gamma(distri_alpha_y+distri_beta_y)*(1/gamma(distri_beta_y))*(1/gamma(distri_alpha_y)),4);
% 
% % exp_meanx=distri_alpha_x/(distri_alpha_x+distri_beta_x);
% % exp_meany=distri_alpha_y/(distri_alpha_y+distri_beta_y);
% %meanx=vpa(int(pdfx*x,x,0,L),4);
% %meany=vpa(int(pdfy*y,y,0,B),4);
% %fplot(pdfx,[0,L]);
% 
% 
% pdf=vpa(pdfx*pdfy,4);
% %exp_mean=exp_meanx*exp_meany;
% %mean=vpa(int(int(x*y*pdf,[0,L]),[0,B]),4);
% %fsurf(pdf,[0 L 0 B]);

 %%
 
count=0;
Sx1=cos(angx1*(pi/180));       
Sy1=cos(angy1*(pi/180));
Sx2=cos(angx2*(pi/180));       
Sy2=cos(angy2*(pi/180));
Sz1=sqrt(1-(Sx1^2 + Sy1^2));
Sz2=sqrt(1-(Sx2^2 + Sy2^2));

for l1=1:1:length(x)
    for l2=1:1:length(x)

a=x1-x(l1,l2);
b=y1-y(l1,l2);
c=H;
R=sqrt((a^2)+(b^2)+(c^2));
psi=acos(H*(R^(-1)));
theta1=acos(((Sx1*a)+(Sy1*b)+(Sz1*c)).*(R^(-1)));
theta2=acos(((Sx2*a)+(Sy2*b)+(Sz2*c)).*(R^(-1)));
% h1=((ga+1)*(1/(R^2))*(0.5*10e-4)*(1/pi)* cos(theta1)*(cos(psi))^(ga))* ((sign(fov(1)/2-abs(theta1))+1)/2);
% h2=((ga+1)*(1/(R^2))*(0.5*10e-4)*(1/pi)* cos(theta2)*(cos(psi))^(ga))* ((sign(fov(2)/2-abs(theta2))+1)/2);
% part1=((sign(h1-th)+1)/2) + ((sign(h2-th)+1)/2) - ((sign(h1-th)+1)/2) * ((sign(h2-th)+1)/2);
h1=((sign(fov(1)/2-abs(theta1))+1)/2);
h2=((sign(fov(2)/2-abs(theta2))+1)/2);
part1=h1+h2-h1*h2;
count= count + part1;
check=check+1;
      if mod(check,10)==0
        check
      end
    end
end

 xx=vpa(subs(count,{angx1, angy1, angx2, angy2},{96.9231  , 96.9231, 73.8462,  101.5385}),5);
 totalre=((cubedim/position_sensi)+1)^2;
 coverage=xx/totalre
 
%%

%%
%  jjj=0;
%  sett=10;
% for l1=50:3:180
%     if (jjj==100)
%         break
%     end
%     for l2=50:3:180
%         if (jjj==100)
%         break
%         end
%         for l3=50:3:180
%             if (jjj==100)
%                      break
%              end
%             for l4=50:3:180
%                 if (jjj==100)
%                       break
%                 end
%                xx=vpa(subs(f,{angx1, angy1, angx2, angy2},{l1  , l2, l3 , l4}),5);
%  
%                  if (xx>sett) 
%                      so1=l1;
%                      so2=l2;
%                      so3=l3;
%                      so4=l4;
%                      sett=xx;
%                  end
%    l4
%             end
%         end
%     end
%     l1
% end
% so1
% so2
% so3
% so4

%  eq1=diff(f,angx1);
% toc   % 38 sec
% eq2=diff(f,angy1);
% toc  %  76 sec
% eq3=diff(f,angx2);
% toc  %   114 sec
% eq4=diff(f,angy2);
% toc  %   152 sec
% jjj=0;
% for l1=50:3:180
%     if (jjj==100)
%         break
%     end
%     for l2=50:3:180
%         if (jjj==100)
%         break
%         end
%         for l3=50:3:180
%             if (jjj==100)
%                      break
%              end
%             for l4=50:3:180
%                 if (jjj==100)
%                       break
%                 end
%                 gg1=vpa(subs(eq1,{angx1, angy1, angx2, angy2},{l1  , l2, l3 , l4}),5);
%                  gg2=vpa(subs(eq2,{angx1, angy1, angx2, angy2},{l1  , l2, l3 , l4}),5);
%                  gg3=vpa(subs(eq3,{angx1, angy1, angx2, angy2},{l1  , l2, l3 , l4}),5);
%                  gg4=vpa(subs(eq4,{angx1, angy1, angx2, angy2},{l1  , l2, l3 , l4}),5);
%                  chh=gg1+gg2+gg3+gg4;
%                  if (chh==0) 
%                      so1=l1;
%                      so2=l2;
%                      so3=l3;
%                      so4=l4;
%                      jjj=100;
%                  end
%    l4
%             end
%         end
%     end
%     l1
% end
% so1
% so2
% so3
% so4
                 
                
% eqns = [eq1 == 0, eq2== 0,eq3 == 0, eq4== 0];
% toc  %   155 sec
% vars=[angx1 angy1 angx2 angy2];

% toc  %   157 sec
%  [solx1, solx2, solx3, solx4] = solve(eqns, vars);
% toc
% solx1=vpa(solx1,3);
% solx2=vpa(solx2,3);
% solx3=vpa(solx3,3);
% solx4=vpa(solx4,3);

%size(solx4)
 
 
% h2=((ga+1)*(1/(R^2))*(0.5*10e-4)*(1/pi)* cos(theta2)*(cos(psi))^(ga));
% h=sqrt(h1^2 + h2^2);
% hnet=vpa(h-threshold,4);
% 
% d_angx1=diff(hnet,angx1);
% d_angx1=vpa(d_angx1,3);
% integrand=vpa(d_angx1*pdf,3);

% dex=0.1;
% dey=0.1;
% check=0;
% exph=0;
% for intx=dex:dex:L
%     check=check+1
%     parsi=0;
%     for inty=dey:dey:B
%         parsi=parsi+(subs(integrand,{x,y},{intx,inty}))*dey;
%     end
%         exph=exph+parsi*dex;
% end
% subs(exph,{angx1, angy1, angx2, angy2},{80,80,100,100})
% exph=vpa(int(int(d_angx1*pdf,x,[0,1]),y,[0,1]),4)

% Ax1=solve(d_angx1==0,angx1);
% Ax1=vpa(Ax1,3);
% exph=vpa(int(int(hnet*pdf,[0,1]),[0,1]),4)

% Ax1=vpa(solve(diff(exph)==0,angx1),4);
% Ax2=vpa(solve(diff(exph)==0,angx2),4);
% Ay1=vpa(solve(diff(exph)==0,angy1),4);
% Ay2=vpa(solve(diff(exph)==0,angy2),4);

toc


