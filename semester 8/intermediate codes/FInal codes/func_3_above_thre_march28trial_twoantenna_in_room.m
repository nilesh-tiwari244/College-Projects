function [enj,new]=func_3_above_thre_march28trial_twoantenna_in_room(x,y,a,b,th,distri)
new=zeros(a,b);
cou=0;
for qq=1:1:a
    for ww=1:1:b
        if(x(qq,ww)> th) || (y(qq,ww)> th)
            new(qq,ww)=1;
            cou=cou+1;
        else
             new(qq,ww)=0;
        end;
    end
end
fin=new.*distri;
enj=sum(sum(fin));
       