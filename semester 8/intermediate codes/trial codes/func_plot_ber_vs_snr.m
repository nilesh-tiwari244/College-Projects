function [] = func_plot_ber_vs_snr(h,low,style)

N=10000;
width=24;
lower_limit=low;
incre=2;
x=randi([0 1],1,N);
mod_data=pskmod(x,2);
BER=zeros(1,width/incre);
new_data=zeros(1,N);

check_plot=0;

for SNR=lower_limit:incre:lower_limit+width
    
check_plot=check_plot+1;
if mod(check_plot,10)==0
        check_plot
    end
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
BER(1,check_plot)=(err1/N);
end
SNR=lower_limit:incre:lower_limit+width;
ylabel('BER');
xlabel('SNR [dB]');
semilogy(SNR(:),BER(1,:),style);
hold on
grid on