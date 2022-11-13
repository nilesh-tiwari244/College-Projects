
function [] = plotbervssnr(h)

N=1000;
width=60;
lower_limit=90;
incre=2.5;


x=randi([0 1],1,N);
mod_data=pskmod(x,2);
BER=zeros(1,width/incre);
new_data=zeros(1,N);

gamma=3;
d=2.5;

check1=0;

for SNR=lower_limit:incre:lower_limit+width
    
check1=check1+1;
err1=0;
en=1/(10^(SNR/10));
 
for tt=1:1:N

noi=normrnd(0,sqrt(en));
recent=mod_data(tt);
new_data(1,:)= noi + h*(recent);

if (new_data(1,tt)>0) && (recent<0)
    err1=err1+1;
elseif (new_data(1,tt)<0) && (recent>0)
       err1=err1+1;
   else 
       err1=err1+0;
   end



end

BER(1,check1)=(err1/N);

end
figure
SNR=lower_limit:incre:lower_limit+width;
ylabel('BER');
xlabel('SNR [dB]');
semilogy(SNR(:),BER(1,:),'-og');
hold on
grid on