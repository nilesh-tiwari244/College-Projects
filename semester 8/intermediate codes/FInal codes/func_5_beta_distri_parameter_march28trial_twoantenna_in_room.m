function [alpha,beta]= func_5_beta_distri_parameter_march28trial_twoantenna_in_room(mode,var)
 alpha=0;
beta=0;
     k=0;   
     b=[10.4/4 8.2/4 5.6/4 4.1/4 3/4 2.3/4 1.8/4 1.3/4 1/4];


     if (var<5.1)
         if (mode<0.41)            
         startwith1=20/var;
         startwith2= startwith1*b(mode*10);
         else
         startwith1=20/var;
         startwith2=startwith1*b(mode*10);
         end
     else
       if (mode<0.41)            
         startwith1=20/var;
         startwith2= startwith1*b(mode*10);
         else
         startwith1=40/var;
         startwith2=startwith1*b(mode*10);
         end
     end
     
%      
% startwith1=1;
% startwith2=1;
         
  for L1=startwith1:0.1:200
     for L2=startwith2:0.1:200
    sensi=0.01;     
 x=0:sensi:1;
 y=zeros(1,length(x));
 index=sensi*2;
 for kk=2:1:length(x)
     y(1,kk)=(betacdf(index,L1,L2)-betacdf((index-sensi),L1,L2));
     index=index+sensi;
 end
 maxy=(max(y));
k=find(y==maxy);
if abs(x(k)-mode)==0
    alpha=L1;
    beta=L2;
    break
end
     end
    if(alpha~=0)
               break
    end
end