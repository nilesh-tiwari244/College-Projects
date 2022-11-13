clc;
clear all;
close all;
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c = averun2;
N=400000;  
x=randi([0 7],1,N);
jj=norm(c);
[m,i] = max(c);
mod_data=pskmod(x,8);
Tx_final = conv(mod_data,c);
pp=length(Tx_final);
limit=30;
BER=zeros(1,limit);
rr=zeros(1,pp);
co=0;
for SNR=0:0.5:limit
  co=co+1;  
    rr=awgn(Tx_final,SNR,'measured');

Rx=rr/jj;
received=zeros(1,length(x));
for q= 1:1:length(x)
    received(q)= Rx(i-1+q);
end
Recovered_bits=pskdemod(received,8);
correct_bits=0;
   diff=Recovered_bits-x;
   for zx=1:1:length(x(1,:))
           if (diff(1,zx)==0)
               correct_bits=correct_bits+1;
           else
               correct_bits=correct_bits+0;
           end
   end
err=((N)-correct_bits);

BER(co)=err/N;
end
figure
SNR=0:0.5:limit;
semilogy(SNR,BER,'-or');
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
% close all;
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c = averun2;
N=400000;  
x=randi([0 3],1,N);
jj=norm(c);
[m,i] = max(c);
mod_data=pskmod(x,4);
Tx_final = conv(mod_data,c);
pp=length(Tx_final);
limit=30;
BER=zeros(1,limit);
rr=zeros(1,pp);
co=0;
for SNR=0:0.5:limit
    co=co+1;
    rr=awgn(Tx_final,SNR,'measured');

Rx=rr/jj;
received=zeros(1,length(x));
for q= 1:1:length(x)
    received(q)= Rx(i-1+q);
end
Recovered_bits=pskdemod(received,4);
% correct_bits=0;
%    diff=Recovered_bits-x;
%    for zx=1:1:length(x(1,:))
%            if (diff(1,zx)==0)
%                correct_bits=correct_bits+1;
%            else
%                correct_bits=correct_bits+0;
%            end
%    end
% err=((N)-correct_bits);
err=symerr(Recovered_bits,x);

BER(co)=err/N;
end
SNR=0:0.5:limit;
semilogy(SNR,BER,'-*b');
hold on;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
%close all;
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c = averun2;
N=400000;
x = randi(2,[1,N])-1; 
jj=norm(c);
                   %Values Given By SIR   
[m,i] = max(c);

for r = 1:1:length(x);         %DEFINING TRANSMITTED SIGNAL USING BPSK MODULATION
if (x(r)==0)
    mod_data(r) = 1;
else
    mod_data(r) = -1;
end
end

Tx_final = conv(mod_data,c);
pp=length(Tx_final);
limit=30;
BER1=zeros(1,limit);
rr=zeros(1,pp);
co=0;
for SNR=0:0.5:limit
    co=co+1;
    rr=awgn(Tx_final,SNR,'measured');

Rx=rr/jj;
y_new=zeros(1,length(x));
for q= 1:1:length(x)
    y_new(q)= Rx(i-1+q);
end
for j = 1:length(y_new)
     if (y_new(j)>0)
         received_bits(j) = 0;
     else
         received_bits(j) = 1;
     end
end
% diff=received_bits-x;
% Err=sum(abs(diff));
err=symerr(received_bits,x);
BER1(co)=err/N;

end
SNR=0:0.5:limit;
semilogy(SNR,BER1,'-+g');
grid;
title('BER.VS.SNR Without OFDM,(for diffrent modulation techniques)');
ylabel('BER');
xlabel('SNR [dB]');
legend('8-Psk','Qpsk','Bpsk');

