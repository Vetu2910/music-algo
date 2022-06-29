%% Introduction

% name:         ARITRA BASU
% stream:       M.tech ECE
% roll :        21ec06008

%% initializations

close all;
clc;
clear all;

j = sqrt(-1);
% segment length
M = 1000;

% segment shift
D = 500;

 
[in,fs] = audioread('taunt.wav');

for i=1:90000
    mix(i)=in(i);
end
% t = linspace(0, 100, 9000);
% fs = 1; %44100
% s1 = sin(3*t)'; %Trying to seperate this and s2
% %[s,fs] = audioread('cymbal_recording_clip.mp3');
% %s_2 = s(:,1);
% s2 = sin(7*t)'; %s_2(3000:3999)';
% mix = s1'+ 0.3*s1'+ 0.7*s2';

mix = mix';

[acf,lags] = autocorr(mix);
Rxx=toeplitz(acf);
%plot(acf)

%% mathworks code for finding the frequency components


% time=mix(:,1);
% signal=mix(:,2);
% h1=spectrum.welch;
% set(h1,'Windowname','Hann');
% Fs=1000;
% set(h1,'OverlapPercent',66.7);
% set(h1,'SegmentLength',2048);
% myPsd=psd(h1,mix(:,2)-mean(mix(:,2)),'Fs',fs)
% semilogx(myPsd.Frequencies,myPsd.Data);xlabel('Frequency(rad/sec)');

%% initializing p and eliminating the necessary eigen values

p = 10;
[coeff1,D] = eig(Rxx);
coeff = fliplr(coeff1);
[row col] = size(coeff);
% so we will take the first 10 eigen values

FV = coeff(:,1:p);     % taken the first p eigen values/vector
noise_subspace = coeff(:,p+1:row);          % creating the noise subspace 
[m n] = size(noise_subspace);
%% creating the complex vector
% jmath = sqrt(-1);
% fq = 1000;
% 
% for i = 1:row
%     s(i)=exp(2*jmath*pi*fq*(i-1));
% end
% s = s';

%% do the fft of the noise vector and then multiply with s vector


noise_fft = fft2(noise_subspace);
noise_fft_mag = abs(noise_fft);

count = 1;
jmath = sqrt(-1);
%fq = 1000;

% problem in s 
w_axis = 0:0.01:pi;
for fq = 1:length(w_axis)
    s = zeros(m,1);
    for i = 1:row
        s(i)=exp(1j*w_axis(fq)*(i-1));  % here is the problem 
    end
% test1 = test1';
% power = ((exp(test1)'.*noise_fft));
 pq = 0;
    for i = 1:n
        pq = pq + abs(s' * noise_subspace(:,i))^2;  % change it 
    end
muzic(count) = 1/pq;
count = count + 1;
% pq = 0;
end
plot(w_axis/pi,10*log10(muzic),'-g');


