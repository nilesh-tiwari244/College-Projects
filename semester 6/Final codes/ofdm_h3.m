clc;
clear all;
close all;

%%  loading CIR and finding the peak point, also calculating the norm value of CIR 
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c= averun2;
% stem(c);   %%%% EXPECTED CHANNEL IMPULSE RESPONSE
jj=norm(c);                    
[m,i] = max(c);
%% Initializing parameters for ofdm part
L1=4096;
L2=(L1)-1;
nor=256;
L3=(nor+1)*2;
Ncp=12;
%% Transmitter

% data generation
Tx_data=randi([0 1],L2,nor);
% Tx_data=randi([0 1],nor,L2);
data_size=nor*L2;
% modulation 
mod_data=pskmod(Tx_data,2);
tx_data1=zeros(L2,L3);
% forming hermitian
for qq=1:1:L2     
    for ww=1:1:nor
        tx_data1(qq,ww+1)=mod_data(qq,ww);
        tx_data1(qq,L3-ww+1)=mod_data(qq,ww)';
    end
end
% Serial to Parallel
s2p=tx_data1.';
% IFFT
am=ifft(s2p);
% convolving channel impulses with input data
vlc5=zeros(length(am(:,1))+length(c)-1,L2);     
for qq=1:1:L2
    vlc5(:,qq)=conv(c,am(:,qq));  
end
% Parallel to series
p2s=vlc5.';
% Cyclic Prefixing
CP_part=p2s(:,end-Ncp+1:end); 
cp=[CP_part p2s];     % appending the cyclic prefix part

%%  Reciever

% Adding Noise using AWGN
SNRstart=0;
SNRincrement=0.5;
SNRend=30;
c=0;
r=zeros(size(SNRstart:SNRincrement:SNRend));
for snr=SNRstart:SNRincrement:SNRend
    c=c+1;
    noisy=awgn(cp,snr,'measured');
    noisy=noisy/jj;            % Normalising
% Removing cyclic prefix part
    cpr=noisy(:,Ncp+1:end); 
    % serial to parallel
     asd2=cpr.';
    % sampling after the peak starts on
     n1=zeros(size(am));
     for qwe=1:1:length(am(1,:))
        for q= 1:1:length(am(:,1))
             n1(q,qwe)= asd2(i-1+q,qwe);
        end
     end

    
   
  %% taking fft and the parallel to serial and then demodulating
% FFT
    amdemod=fft(n1);
% Parallel to serial
    rserial=amdemod.';
% QAM demodulation 
    Recovered_bits=pskdemod(rserial,2);
    
    
%% Calculating the Bit Error Rate
Recovered_bits=Recovered_bits(:,2:nor+1);
  err=symerr(Recovered_bits,Tx_data);

BER(c)=err/data_size;

end
snr=SNRstart:SNRincrement:SNRend;
%% Plotting BER vs SNR
semilogy(snr,BER,'-+g');
hold on;
grid;
title('BER.VS.SNR, with OFDM(for diffrent modulation techniques)');
ylabel('BER');
xlabel('SNR [dB]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;

%%  loading CIR and finding the peak point, also calculating the norm value of CIR 
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c= averun2;
% stem(c);   %%%% EXPECTED CHANNEL IMPULSE RESPONSE
jj=norm(c);                    
[m,i] = max(c);
%% Initializing parameters for ofdm part
L1=4096;
L2=(L1)-1;
nor=256;
L3=(nor+1)*2;
Ncp=12;
%% Transmitter

% data generation
Tx_data=randi([0 3],L2,nor);
% Tx_data=randi([0 1],nor,L2);
data_size=nor*L2;
% modulation 
mod_data=pskmod(Tx_data,4);
tx_data1=zeros(L2,L3);
% forming hermitian
for qq=1:1:L2     
    for ww=1:1:nor
        tx_data1(qq,ww+1)=mod_data(qq,ww);
        tx_data1(qq,L3-ww+1)=mod_data(qq,ww)';
    end
end
% Serial to Parallel
s2p=tx_data1.';
% IFFT
am=ifft(s2p);
% convolving channel impulses with input data
vlc5=zeros(length(am(:,1))+length(c)-1,L2);     
for qq=1:1:L2
    vlc5(:,qq)=conv(c,am(:,qq));  
end
% Parallel to series
p2s=vlc5.';
% Cyclic Prefixing
CP_part=p2s(:,end-Ncp+1:end); 
cp=[CP_part p2s];     % appending the cyclic prefix part

%%  Reciever

% Adding Noise using AWGN
SNRstart=0;
SNRincrement=0.5;
SNRend=30;
c=0;
r=zeros(size(SNRstart:SNRincrement:SNRend));
for snr=SNRstart:SNRincrement:SNRend
    c=c+1;
    noisy=awgn(cp,snr,'measured');
    noisy=noisy/jj;            % Normalising
% Removing cyclic prefix part
    cpr=noisy(:,Ncp+1:end); 
    % serial to parallel
     asd2=cpr.';
    % sampling after the peak starts on
     n1=zeros(size(am));
     for qwe=1:1:length(am(1,:))
        for q= 1:1:length(am(:,1))
             n1(q,qwe)= asd2(i-1+q,qwe);
        end
     end

    
   
  %% taking fft and the parallel to serial and then demodulating
% FFT
    amdemod=fft(n1);
% Parallel to serial
    rserial=amdemod.';
% QAM demodulation 
    Recovered_bits=pskdemod(rserial,4);
    
    
%% Calculating the Bit Error Rate
Recovered_bits=Recovered_bits(:,2:nor+1);
  err=symerr(Recovered_bits,Tx_data);

BER(c)=err/data_size;

end
snr=SNRstart:SNRincrement:SNRend;
%% Plotting BER vs SNR
semilogy(snr,BER,'-*b');
hold on;
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;

%%  loading CIR and finding the peak point, also calculating the norm value of CIR 
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c= averun2;
% stem(c);   %%%% EXPECTED CHANNEL IMPULSE RESPONSE
jj=norm(c);                    
[m,i] = max(c);
%% Initializing parameters for ofdm part
L1=4096;
L2=(L1)-1;
nor=256;
L3=(nor+1)*2;
Ncp=12;
%% Transmitter

% data generation
Tx_data=randi([0 7],L2,nor);
% Tx_data=randi([0 1],nor,L2);
data_size=nor*L2;
% modulation 
mod_data=pskmod(Tx_data,8);
tx_data1=zeros(L2,L3);
% forming hermitian
for qq=1:1:L2     
    for ww=1:1:nor
        tx_data1(qq,ww+1)=mod_data(qq,ww);
        tx_data1(qq,L3-ww+1)=mod_data(qq,ww)';
    end
end
% Serial to Parallel
s2p=tx_data1.';
% IFFT
am=ifft(s2p);
% convolving channel impulses with input data
vlc5=zeros(length(am(:,1))+length(c)-1,L2);     
for qq=1:1:L2
    vlc5(:,qq)=conv(c,am(:,qq));  
end
% Parallel to series
p2s=vlc5.';
% Cyclic Prefixing
CP_part=p2s(:,end-Ncp+1:end); 
cp=[CP_part p2s];     % appending the cyclic prefix part

%%  Reciever

% Adding Noise using AWGN
SNRstart=0;
SNRincrement=0.5;
SNRend=30;
c=0;
r=zeros(size(SNRstart:SNRincrement:SNRend));
for snr=SNRstart:SNRincrement:SNRend
    c=c+1;
    noisy=awgn(cp,snr,'measured');
    noisy=noisy/jj;            % Normalising
% Removing cyclic prefix part
    cpr=noisy(:,Ncp+1:end); 
    % serial to parallel
     asd2=cpr.';
    % sampling after the peak starts on
     n1=zeros(size(am));
     for qwe=1:1:length(am(1,:))
        for q= 1:1:length(am(:,1))
             n1(q,qwe)= asd2(i-1+q,qwe);
        end
     end

    
   
  %% taking fft and the parallel to serial and then demodulating
% FFT
    amdemod=fft(n1);
% Parallel to serial
    rserial=amdemod.';
% QAM demodulation 
    Recovered_bits=pskdemod(rserial,8);
    
    
%% Calculating the Bit Error Rate
Recovered_bits=Recovered_bits(:,2:nor+1);
  err=symerr(Recovered_bits,Tx_data);

BER(c)=err/data_size;

end
snr=SNRstart:SNRincrement:SNRend;
%% Plotting BER vs SNR
semilogy(snr,BER,'-or');
hold on;
grid;
legend('Bpsk','Qpsk','8-Psk');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

