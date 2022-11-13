l=3;
limit=l^(-1.5);
syms y;
part1=sqrt((y^(-0.3333))-{l^2});
fd=(y^(-4/3)) * (1/part1) *  (part1/4) * (1/exp((part1^2)/8)) ;
k1=double(int(fd,y,0,limit));
cd=1/k1;
