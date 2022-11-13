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
L=2048;
L1=L/2;
L2=(L1)-1;
nor=96;
Ncp=12;
%% Transmitter

% data generation
Tx_data=randi([0 7],nor,L2);
% Tx_data=randi([0 1],nor,L2);
data_size=nor*L2;
% modulation 
 mod_data=pskmod(Tx_data,8);
tx_data1=zeros(nor,L);
% forming hermitian
for qq=1:1:nor     
    for ww=1:1:L2
        tx_data1(qq,ww+1)=mod_data(qq,ww);
        tx_data1(qq,L-ww+1)=mod_data(qq,ww)';
    end
end
% Serial to Parallel
s2p=tx_data1.';
% IFFT
am=ifft(s2p);
% convolving channel impulses with input data
vlc5=zeros(length(am(:,1))+length(c)-1,nor);     
for qq=1:1:nor
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
SNRincrement=1;
SNRend=15;
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
Recovered_bits=Recovered_bits(:,2:L2+1);
  err=symerr(Recovered_bits,Tx_data);

BER(c)=err/data_size;

end
snr=SNRstart:SNRincrement:SNRend;
%% Plotting BER vs SNR
semilogy(snr,BER,'-or');
hold on;
grid;
title('BER.VS.SNR, with OFDM (for diffrent modulation techniques)');
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
L=2048;
L1=L/2;
L2=(L1)-1;
nor=96;
Ncp=12;
%% Transmitter

% data generation
Tx_data=randi([0 3],nor,L2);
data_size=nor*L2;
% modulation 
 mod_data=pskmod(Tx_data,4);
tx_data1=zeros(nor,L);
% forming hermitian
for qq=1:1:nor     
    for ww=1:1:L2
        tx_data1(qq,ww+1)=mod_data(qq,ww);
        tx_data1(qq,L-ww+1)=mod_data(qq,ww)';
    end
end
% Serial to Parallel
s2p=tx_data1.';
% IFFT
am=ifft(s2p);
% convolving channel impulses with input data
vlc5=zeros(length(am(:,1))+length(c)-1,nor);     
for qq=1:1:nor
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
SNRincrement=1;
SNRend=15;
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
Recovered_bits=Recovered_bits(:,2:L2+1);
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
s

%%  loading CIR and finding the peak point, also calculating the norm value of CIR 
addpath('C:\Users\User\Documents\6thsembtp\study material\by sir\CIRs\Scenario 3\Home\D3');
load('Run1.mat');
c= averun2;
% stem(c);   %%%% EXPECTED CHANNEL IMPULSE RESPONSE
jj=norm(c);                    
[m,i] = max(c);
%% Initializing parameters for ofdm part
L=2048;
L1=L/2;
L2=(L1)-1;
nor=96;
Ncp=12;
%% Transmitter

% data generation
Tx_data=randi([0 1],nor,L2);
% Tx_data=randi([0 1],nor,L2);
data_size=nor*L2;
% QAM modulation 
%mod_data=qammod(Tx_data,16);
 mod_data=pskmod(Tx_data,2);
tx_data1=zeros(nor,L);
% forming hermitian
for qq=1:1:nor     
    for ww=1:1:L2
        tx_data1(qq,ww+1)=mod_data(qq,ww);
        tx_data1(qq,L-ww+1)=mod_data(qq,ww)';
    end
end
% Serial to Parallel
s2p=tx_data1.';
% IFFT
am=ifft(s2p);
% convolving channel impulses with input data
vlc5=zeros(length(am(:,1))+length(c)-1,nor);     
for qq=1:1:nor
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
SNRincrement=1;
SNRend=15;
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
Recovered_bits=Recovered_bits(:,2:L2+1);
  err=symerr(Recovered_bits,Tx_data);

BER(c)=err/data_size;

end
snr=SNRstart:SNRincrement:SNRend;
%% Plotting BER vs SNR
semilogy(snr,BER,'-+g');
hold on;
grid;

legend('8-Psk','Qpsk','Bpsk');