function [pp]=func_6_generate_random(distri,first,last)
len=length(distri);
 x=rand*(len);
%x=rand*(1000000);
range=zeros(2,len);
count=0;
for l1=1:1:len
    range(1,l1)=count;
    range(2,l1)=count+len*distri(l1);
     % range(2,l1)=count+1000000*distri(l1);
    count=range(2,l1);
end
for l1=1:1:len
%     if((x<range(2,l1)) && (x>range(1,l1)))
 if (x<range(2,l1))
        pp=first + (l1-1)*(((last-first))/(len-1));
        break
    end
end
