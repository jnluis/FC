%% Problema 6.1 - Propriedades da transformada de Fourier discreta
clc;
clear all;
close all;

h=0.1;
N=2^10; % 1024. fft e ifft funcionam melhor com potências de 2
t=0:h:(N-1)*h;  

homg=2*pi/(N*h); % = Delta omega
% Com shift, os valores omega são:
omgmax=(N/2-1)*homg;
omgmin=(-N/2)*homg;
omg=omgmin:homg:omgmax;

nyquist=pi/h % Frequencia de Nyquist
    
% Equação

y = exp(-10*t*1i) + exp(1i*20*t);
z= fftshift(fft(y)); 
dens_esp=(h*abs(z)).^2;

subplot(2,1,1)
plot(omg, dens_esp)
title("\it{y}=exp(-i10t) + exp(i20t)")

subplot(2,1,2)
plot(t, abs(y))