%% Introduction

% name:         ARITRA BASU
% stream:       M.tech ECE
% roll :        21ec06008

%% Initializations

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
length_audio = length(mix');

% no of segments 

L = (length_audio/D)-1;

%% segmenting the signals

bw = bartlett(M);
% plot(bw);
U = sum(abs(bw).^2)*(1/M);
% pxx = zeros(1025,179);
start = 1;
count = 1;

for j = 1:L
    for i = start:(start+M-1)
        signal_in(count)=(mix(i)*bw(count)*exp(-2*j*pi*fs*count));   % multiplying by the window function
        signal_in1(count) = mix(i)*bw(count);
        count = count + 1;
    end

%signal_out(j) = norm(fftshift(fft(signal_in)));
%pxx(:,j) = periodogram(signal_in)/U;
signal_out1(j,:)=fft(signal_in1).^2/(M*U);
signal_out(j,:) = (signal_in).^2/(M*U);
start = start + D;
count = 1 ;
end

%% finding the welsh PSD estimate 

% signal_out = signal_out/(M*U);
pxx_final = sum(abs(signal_out1))/L;
figure;
plot(abs(pxx_final));
figure;
plot(mix);
