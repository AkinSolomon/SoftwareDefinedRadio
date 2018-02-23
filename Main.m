%% Software Template
clc
clear all

EbperN0dB =10; % SNR
EbperN0 = 10^(EbperN0dB/10);

alpha = 0.5; % SRRC rolloff param
Fs=15000;  %Sampling Frequency
D = 5; % truncation to [-DT,DT]
Ts=1/Fs; % Sampling Rate
L = floor(0.0025/Ts); % oversampling factor
H=2;                 % Gain and Phase Offset
pilots=[1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1]; % pilot type
%%%start = 40;% where to place markers into noise
M=2; %constellation size; "2" for BPSK


% Unknown Values, Guess to make program run
fc=3000;  %Carrier Frequency
fdel=0; %Frequency offset, fdel=phase_error in radians/ 2piT

pb=900;%0.5*(Fs/L)*(1+alpha); %1500;    %Passband edge
ps=2*fc - pb; %4500;    %Stopband edge
%%%phi=0;       %phase offset
% %% Table 1.1, step (a)
% % Create an audiorecorder object with 8000 sps and
% % a single channel (mono); view its properties:
% recObj = audiorecorder(Fs, 16, 1);
% get(recObj);
% 
% % Collect a sample of your speech with a microphone;
% % and, plot the signal data:
% 
% % Record your voice for 2 seconds. Use display and 
% 
% % pause as aids to control the start of the recording.
% Trec = 3; %2 second record time
% disp('Press Enter to start recording.')
% pause;%wait for keystroke
% disp('Recording.')
% recordblocking(recObj, Trec);
% disp('End of Recording.')
% 
% % Store data in double-precision array.
% a1 = getaudiodata(recObj);
% 
% %Convert 8bit rep to symbols
% % a1(a1==0)=1;
% 
% % Plot the waveform.
% t = (0:length(a1)-1)/Fs; %sample times (sec)
% plot(t,a1);
% 
% N=length(a1);



%% This section takes in the input data source and converts
% it into a string of bits.
prompt='Enter Data Source (r: random bits/ t:text input):';

x=input(prompt,'s');
M=input('Enter Modulation Scheme (2:BPSK, 4:QPSK):');
varw =(EbperN0*log2(M))^(-1)/2; %per channel sample var

switch x
    case 'r'
        N = 235; %number of symbols
        if(M==2)
            a = (2*round(rand(1,N))-1); %noiseless I
            Nbits=length(a)*2;
            
        elseif(M==4)
            a = (2*round(rand(1,N))-1) ... %noiseless I
            + 1j*(2*round(rand(1,N))-1); %and add Q
            a = a/sqrt(2); %unit length QPSK symbols
        end
    case 't'
        text=input('Enter input string:','s');
        [a, bits2]= char2psk(text,M);
        N=length(bits2);
  end

% Channel Impairments: Visualize Transmitted bits
i=1;
% 
a1=[a pilots];% postpend pilots
a=a1';
while (mod(length(a),8)~=0)
a=[a zeros(1,i)];
end

y=a;

start=2;
v = sqrt(varw/2)*randn(1,3*length(a1));
indices = start: start+length(a1)-1;
v(indices)=v(indices)+ H*a1;%add scaled baseband message into noise
y=v;

figure;
plot(real(y),imag(y),'x','MarkerSize',8,'LineWidth',2);
title('Scatter Plot','fontsize',13)
xlabel('In-phase','fontsize',13);
ylabel('Quadrature','fontsize',13);
axis square% set axes extents
grid % draw grid lines

% Pulse Shaping: This section maps the transmitted symbols unto an analog waveform / Modulation and Transmission
 y1=upsample(y,L); 
  g1=srrc(D,alpha,L);
 
 y_up=conv(y1,g1,'same'); % Output of Analog Waveform
 
 
% Visualize Eye Diagram
 eyediagram(y_up,L,N,'complex')
 

 
 
 
%% Modulation
t2=0:1/Fs:(length(y_up)-1)/Fs;
s=exp(1j*2*pi*(fc)*t2);

m2t=y_up.*s;

h=firlpf(48,pb,ps,Fs);

plottf(real(m2t),Ts,'f')


 %% Demodulation
 t2=0:1/Fs:(length(y_up)-1)/Fs;
 s=2*exp(-1j*2*pi*(fc+fdel)*t2);
 %%%omit fdel at receiver, or set to zero (useful in simulation only)
 
 
%  m3t=real(m2t); 
 m2t=conv(m2t.*s,h, 'same');
 
 m4t=conv(m2t,g1);%,'same'); %Match with reciever filter

 figure;plottf(real(m4t),Ts);%,'f')
%  eyediagram(m4t,L,N,'complex')
 
figure;subplot(211);plot(real(m4t));subplot(212);plot(imag(m4t))
figure;plot(abs(m4t))
 %% Packet Detection
 
 len_pak=length(pilots)+50*8/log2(M); %50 ASCII characters payload
 m4t=packetdetect(m4t,L,len_pak);
 %%%13 is hard coded adn depend on pilots; replace with length(pilots)
 %%%60 is hard coded and matches teh number of characters in the phoen app
 
 
 
% plot(real(m2t),imag(m2t));
% L2=10;
% hold on;
% plot(real(m2t(1:L2:end)),imag(m2t(1:L2:end)),'ro');
% 
% hold off; title('IQ Plot')
% xlabel('In-phase');
% ylabel('Quadrature');



%% Symbol Timing
% CMA
% for p=0:L-1
% z(p+1)=var(abs(m4t(p + 1:L:end)));
% end
% 
% [min_val,ind]=min(z);
% p_est=ind-1;
% 
% plot(z)

v=m4t;
% Estimate Symbol timing
cost = zeros(1,L);
%v is the recieved symbols
for p = 0:L-1
    gamma_vector = ones(1,length(v(1+p:L:length(v))));
    gamma = mean(abs(v(1+p:L:length(v))));
    gamma_vector = gamma*gamma_vector;
    cost(p+1) = sum((abs(v(1+p:L:length(v)))-gamma_vector).^2);
end

[~,offset] = min(cost);
offset = offset -1;


y1=m4t(1:L:end);
y=m4t(offset:L:end); %downsampling

L2=10;
figure;
plot(real(y(1:L2:end)), ...
imag(y(1:L2:end)),'bo');

pause(1)

hold on
plot(real(y1(1:L2:end)), ...
imag(y1(1:L2:end)),'ro');
grid on; axis square;
title('IQ samples')

plot(real(a(1:L2:end)), ...
imag(a(1:L2:end)),'gd');
grid on; axis square;
title('IQ samples')

% ylim([-1 1])
% xlim([-3 3])

hold off

legend('After symbol timing','Before symbol timing','Original Message')

%% Frame recovery
xc=conv(y,fliplr(pilots),'same'); %cross-correlation with pilots

[value, indx]=max(abs(xc));
peak=xc(indx);

H_est=peak/length(pilots);

y1=y./H_est;

% correlation output
figure;plot(abs(xc));
xlabel('lag number');title('crosscorrelation')


% Cut frame from noisy data

%%%insert here
Npak=13+50*8/log2(M);
start = indx-length(pilots)+1;

symbolsRX=y1(start:1:start+Npak-1);

%%% use name "symoblsRX or better name of your choicee to match below  

figure;

plot(real(a(1:L2:end)), ...
imag(a(1:L2:end)),'gd');
grid on; axis square;
title('IQ samples')

hold on
pause(2)% view at two frames per second


y_tmp=y1(symbolsRX);

plot(real(y_tmp(1:L2:end)), ...
imag(y_tmp(1:L2:end)),'ro');
grid on; axis square;
title('IQ samples')


ylim([-1.5 1.5])


hold off
legend( 'Original Message','After frame timing')



%% Frequency Recovery

r=angle(pilots./y1(indx-6:indx+6));
%%%will be indices 1:13

r=unwrap(r);

p=polyfit(pac_len(1:13),r,1);
%l=r-polyval(p,y1(indx-6:indx+6));

f_d=p(1)/(2*pi*L/Fs);
phi=p(2);

%plot(1:length(pilots),l)
% plots
hold on
plot(1:length(pilots),r, 'rd')
grid on
hold off

phase_err=exp(1j*(2*pi*f_d*L/Fs*pac_len+ phi));
%%% yes, need to be consistent with indexind here, to continue straight line
a_rec=y1(pac_len).*phase_err;

plot(real(a_rec),imag(a_rec),'ro')
hold on
plot(real(y_tmp),imag(y_tmp),'bd')
title('IQ Plot before and after freq recovery')

legend('After frequency recovery','Before Frequency Recovery')


%% Visualization
phasediff=angle(a./a_rec);
figure; plot(phasediff);
axis square
title('Phase variation across the samples')
 
 %% Converts Recieved Signal back into bits.
 
[m,Txbits]=psk2char(a,2);
[m,Rxbits]=psk2char(a_rec,2);


BER=length(find(Rxbits-Txbits))/N;

 
 
 