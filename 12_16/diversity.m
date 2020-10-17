close all;
clear all;

%% params
M=2; % M-ary modulation - M different symbols
n=10000; % n data symbols to be trasnmitted
T_s=1e-5; % 10 us (10^-5 second) symbol duration, 100 KHz signal bandwidth
f_s = 1/T_s; % frequency rate
T_sample=1e-6; %Sampling duration = 0.25us

ricChan = comm.RicianChannel('SampleRate',1e5, 'MaximumDopplerShift',10, 'KFactor',2);
rayChan = comm.RayleighChannel('SampleRate',1e5,'MaximumDopplerShift',10);
%rayChan.DopplerSpectrum = doppler('Gaussian',1);
testChan = rayChan;

%% plot rayleigh-channel
figure(1);
ch = rayChan(ones(1e5,1));
db_ch = 10*log10(abs(ch));
plot(0:T_s:1-T_s,db_ch)
title('Rayleigh Channel');
ylabel('Received Power (dB)');
xlabel('Time (s)');


data=randi(M,1,n)-1; %%generate n random data symbols

%I and Q components of the intended data signal
tx = repelem(data,round(T_s/T_sample))';

%use DBPSK to encode data
mod = comm.DBPSKModulator;
%use BPSK to decode data
demod = comm.BPSKDemodulator;

% 使傳輸的訊號透過rayley的channel
fadedSig = testChan(mod(tx));

%套用特定的SNR值，把noise加上去
awgnChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (SNR)');

% 計算bit error rate
errorCalc = comm.ErrorRate;

SNR = 0:2:20; % Range of SNR values, in dB.
numSNR = length(SNR);
berVec = zeros(3, numSNR); % Preallocate a vector for BER results 

% 透過設定target的SNR來方便操作
for n = 1:numSNR
   awgnChan.SNR = SNR(n);
   rxSig = awgnChan(fadedSig); % Add Gaussian noise
   %% Manually convert DBPSK to BPSK
   pf = [1 ;rxSig(1:end-1)];
   rxSig = rxSig./pf;
   
   rx = demod(rxSig);  % Demodulate
   reset(errorCalc)
   berVec(:,n) = errorCalc(tx,rx); % Compute error rate.
end
%bit error rate
BER = berVec(1,:);

BERtheory = berawgn(SNR,'dpsk',M,1);

% 將有faded和沒有faded的做一個比較
figure(2);
hold on;
plot(SNR,db(BERtheory),'b-',SNR,db(BER),'r*');
legend('Theoretical AWGN BER','Ray BER');
xlabel('SNR (dB)'); ylabel('BER (db)');
ylim([-50 0]);
title('DBPSK over Rayleigh Fading Channel');

%{
%% Selection Combining (SC)

fadedSig2 = testChan(mod(tx));  %simulate the signal from second antenna
rx_SC = zeros(size(tx));
berVec_SC = zeros(3, numSNR);
for n = 1:numSNR
   awgnChan.SNR = SNR(n);
   rxSig = awgnChan(fadedSig); % Add Gaussian noise
   rxSig2 = awgnChan(fadedSig2);
   %% Manually convert DBPSK to BPSK
   idx = abs(rxSig2)>abs(rxSig);
   
   pf = [1 ;rxSig(1:end-1)];
   rxSig = rxSig./pf;
   pf = [1 ;rxSig2(1:end-1)];
   rxSig2 = rxSig2./pf;

   rx_SC = rxSig;
   rx_SC(idx) = rxSig2(idx);
   rx_SC = demod(rx_SC);
   %rx_SC = rx*abs(rxSig)*abs(rxSig)+;
   
   reset(errorCalc);
   berVec_SC(:,n) = errorCalc(tx,rx_SC); % Compute error rate.
end

BER_SC = berVec_SC(1,:);
%}

%% implement 4 branches Selection Combining (SC)

fadedSig2 = testChan(mod(tx));  %simulate the signal from second antenna
fadedSig3 = testChan(mod(tx));  %simulate the signal from third antenna
fadedSig4 = testChan(mod(tx));  %simulate the signal from forth antenna

rx_4SC = zeros(size(tx));
berVec_4SC = zeros(3, numSNR);
for n = 1:numSNR
    
   awgnChan.SNR = SNR(n);
   rxSig = awgnChan(fadedSig); % Add Gaussian noise
   rxSig2 = awgnChan(fadedSig2);
   rxSig3 = awgnChan(fadedSig3);
   rxSig4 = awgnChan(fadedSig4);
   
   %% Manually convert DBPSK to BPSK
   rx_Sig_all = [rxSig'; rxSig2'; rxSig3'; rxSig4'];
   [value, idx] = max(abs(rx_Sig_all));
   
   pf = [1 ;rxSig(1:end-1)];
   rxSig = rxSig./pf;
   pf = [1 ;rxSig2(1:end-1)];
   rxSig2 = rxSig2./pf;
   pf = [1 ;rxSig3(1:end-1)];
   rxSig3 = rxSig3./pf;
   pf = [1 ;rxSig4(1:end-1)];
   rxSig4 = rxSig4./pf;

   rx_Sig_all = [rxSig'; rxSig2'; rxSig3'; rxSig4'];
   for i = 1:size(rx_Sig_all,2)
       rx_4SC(i) = rx_Sig_all(idx(i), i);
   end
   rx_4SC = demod(rx_4SC);
   %rx_SC = rx*abs(rxSig)*abs(rxSig)+;
   
   reset(errorCalc);
   berVec_4SC(:,n) = errorCalc(tx,rx_4SC); % Compute error rate.
end

BER_4SC = berVec_4SC(1,:);

%% implement 2 branches Maximum Ratio Combining (MRC) 

rx_2MRC = zeros(size(tx));
berVec_2MRC = zeros(3, numSNR);
for n = 1:numSNR
    
   awgnChan.SNR = SNR(n);
   rxSig = awgnChan(fadedSig); % Add Gaussian noise
   rxSig2 = awgnChan(fadedSig2);
   
   %% Manually convert DBPSK to BPSK
   sig1_power = abs(rxSig);
   sig2_power = abs(rxSig2);
   sig1_ratio = sig1_power ./ (sig1_power + sig2_power);
   sig2_ratio = sig2_power ./ (sig1_power + sig2_power);
   
   pf = [1 ;rxSig(1:end-1)];
   rxSig = rxSig./pf;
   pf = [1 ;rxSig2(1:end-1)];
   rxSig2 = rxSig2./pf;

   for i = 1:size(rx_Sig_all,2)
       rx_2MRC(i) = sig1_ratio(i) * rxSig(i) + sig2_ratio(i) * rxSig2(i);
   end
   rx_2MRC = demod(rx_2MRC);
   %rx_SC = rx*abs(rxSig)*abs(rxSig)+;
   
   reset(errorCalc);
   berVec_2MRC(:,n) = errorCalc(tx,rx_2MRC); % Compute error rate.
end

BER_2MRC = berVec_2MRC(1,:);

%% implement 4 branches Maximum Ratio Combining (MRC) 

rx_4MRC = zeros(size(tx));
berVec_4MRC = zeros(3, numSNR);
for n = 1:numSNR
    
   awgnChan.SNR = SNR(n);
   rxSig = awgnChan(fadedSig); % Add Gaussian noise
   rxSig2 = awgnChan(fadedSig2);
   rxSig3 = awgnChan(fadedSig3);
   rxSig4 = awgnChan(fadedSig4);
   
   %% Manually convert DBPSK to BPSK
   sig1_power = abs(rxSig);
   sig2_power = abs(rxSig2);
   sig3_power = abs(rxSig3);
   sig4_power = abs(rxSig4);
   
   sig1_ratio = sig1_power ./ (sig1_power + sig2_power +  sig3_power +  sig4_power);
   sig2_ratio = sig2_power ./ (sig1_power + sig2_power +  sig3_power +  sig4_power);
   sig3_ratio = sig3_power ./ (sig1_power + sig2_power +  sig3_power +  sig4_power);
   sig4_ratio = sig4_power ./ (sig1_power + sig2_power +  sig3_power +  sig4_power);
   
   pf = [1 ;rxSig(1:end-1)];
   rxSig = rxSig./pf;
   pf = [1 ;rxSig2(1:end-1)];
   rxSig2 = rxSig2./pf;
   pf = [1 ;rxSig3(1:end-1)];
   rxSig3 = rxSig3./pf;
   pf = [1 ;rxSig4(1:end-1)];
   rxSig4 = rxSig4./pf;

   for i = 1:size(rx_Sig_all,2)
       rx_4MRC(i) = sig1_ratio(i) * rxSig(i) + sig2_ratio(i) * rxSig2(i) + sig3_ratio(i) * rxSig3(i) + sig4_ratio(i) * rxSig4(i);
   end
   rx_4MRC = demod(rx_4MRC);
   %rx_SC = rx*abs(rxSig)*abs(rxSig)+;
   
   reset(errorCalc);
   berVec_4MRC(:,n) = errorCalc(tx,rx_4MRC); % Compute error rate.
end

BER_4MRC = berVec_4MRC(1,:);


figure(3);
hold on;
plot(SNR, db(BER), 'r*', SNR, db(BER_4SC), 'bo', SNR, db(BER_2MRC), 'k.', SNR, db(BER_4MRC), 'm+');
legend('No SC BER', '4-SC BER', '2-MRC BER', '4-MRC BER');
xlabel('SNR (dB)'); ylabel('BER (dB)');
title('DBPSK over Rayleigh Fading Channel');