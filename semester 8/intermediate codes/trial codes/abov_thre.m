function [enj,new,neth,heqv]=abov_thre(x,y,a,b,th,distri)
new=zeros(a,b);
heqv=new;
cou=0;
for qq=1:1:a
    for ww=1:1:b
        heqv(qq,ww)=sqrt((x(qq,ww)^2)+(y(qq,ww)^2));
        if(x(qq,ww)> th) || (y(qq,ww)> th)
            new(qq,ww)=1;
            cou=cou+1;
        else
             new(qq,ww)=0;
        end;
    end
end
neth=sum(sum(heqv.*distri));
fin=new.*distri;
enj=sum(sum(fin));
        