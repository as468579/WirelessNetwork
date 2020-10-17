%%BPSK Narrowband Interference
%This code will generate a QPSK waveform modulated with random data
%Then compare the frequency domain spectrum with and without square root
%raised cosine pulse shaping.

clear;
M=2 % M-ary modulation - M different symbols
n=1000; % n data symbols to be trasnmitted
data=randi(M,1,n); %%generate n random data symbols
f_c=100e3; % 100 KHz carrier frequency
T_s=1e-4; % 100 us (10^-4 second) symbol duration, 10 KHz signal bandwidth
T_sample=0.25e-6; %Sampling duration = 0.25us
inteference_att=3; %attenuation of the interference in dB

T_s_i=1e-2; % symbol duration of the interfering signal
data_i=randi(M,1,n/round(T_s_i/T_s)); %n data symbols in the interfering signal
f_c_i=103e3; % carrier frequency of the interfering signal
phase_i=rand(1)*2*pi;

%%I and Q symbol table (BPSK)
I_table=[1 -1];
Q_table=[0 0];


N_spreading=8; % number of chips per symbol
data_s=repelem(data, N_spreading)-1; %repeat each element for N_spreading times
key = randi(2, size(data_s))-1;

%%!!! implement spreading here!!! (replace the following line)
data_spreaded=data_s; 
data_spreaded = xor(data_s, key);%%%%%%%%%%%%%

data_spreaded = data_spreaded + 1;

%I and Q components of the intended data signal
I=I_table(data_spreaded);
Q=Q_table(data_spreaded);
%I and Q components of the intended interfering signal
I_i=I_table(data_i);
Q_i=Q_table(data_i);



%resample I and Q at sampling frequency
I_m=repelem(I, round(T_s/N_spreading/T_sample));
Q_m=repelem(Q, round(T_s/N_spreading/T_sample));

I_i_m=repelem(I_i, round(T_s_i/T_sample));
Q_i_m=repelem(Q_i, round(T_s_i/T_sample));

t=[0:T_sample:n*T_s-T_sample];


%x=I_m.*cos(2*pi*f_c*t)+Q_m.*sin(2*pi*f_c*t);
%x_i=I_i_m.*cos(2*pi*f_c_i*t+phase_i)+Q_i_m.*sin(2*pi*f_c_i*t+phase_i);

x=((I_m+j*Q_m).*exp(j*2*pi*f_c*t));
x_i=((I_i_m+j*Q_i_m).*exp(j*2*pi*f_c_i*t));
x_i=x_i*10^(inteference_att/10);

figure;


f=[0:1/(n*T_s):1/T_sample-1]-.5*1/T_sample;

f_x_db=10*log10(abs(fftshift(fft(x))).^2);
f_x_i_db=10*log10(abs(fftshift(fft(x_i))).^2);
max_f_x_db=max(f_x_db);
plot(f,f_x_db-max_f_x_db);
hold on;
plot(f,f_x_i_db-max_f_x_db);
axis([f_c-2*(N_spreading/T_s) f_c+2*(N_spreading/T_s) -100 0]);
ylabel('Normalized power (dB)');
xlabel('Frequency (Hz)');
legend('Intented Data Signal','Narrowband Interfering Signal');

%add some noise
noise=(randn(1,length(t))+randn(1,length(t))*j)*0.01;
%received signal = transmitted signal + noise
r=x+noise;
%downconvert to baseband
r=r.*exp(-j*2*pi*f_c*t);

%same for the received signal when considering interference
r_i=x+x_i+noise;
r_i=r_i.*exp(-j*2*pi*f_c*t);


cps=round(T_s/T_sample/N_spreading); %no. of samples per chip



%averaging over each chip duration
r_chip=arrayfun(@(i) mean(r_i(i:i+cps-1)),1:cps:length(r_i)-cps+1)';



%decision device
r_chip(find(real(r_chip)>=0))=1;
r_chip(find(real(r_chip)<0))=2;
%%%%%%%%%%%
r_chip = r_chip-1;
r_chip = xor(r_chip, key');
r_chip = r_chip +1;
%%%%%%%%%%%
r_symbol=r_chip;
%averaging over each symbol duration
r_symbol=arrayfun(@(i) mean(r_symbol(i:i+N_spreading-1)),1:N_spreading:length(r_symbol)-N_spreading+1)';
%decision device
r_symbol(find(real(r_symbol)>=1.5))=2;
r_symbol(find(real(r_symbol)<1.5))=1;

%r_symbol=repelem(r_symbol, N_spreading)-1; %repeat each element for N_spreading times

% r_symbol = xor(r_symbol, key);
% r_symbol = r_symbol +1;

sum(abs(data-r_symbol.'))/n %error rate









