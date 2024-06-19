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

%y = sin(10*t) + sin(10.05*t) 

y = sin(10*t) + sin(10.05*t);
z= fftshift(fft(y)); 
dens_esp=(h*abs(z)).^2;

figure(1)
subplot(2,1,1)
plot(omg, dens_esp)
title("\it{y}=sin(10t) + sin(10.05\it{t})")

subplot(2,1,2)
plot(t, y)

%% Corrigido o valor de N

h=0.1;
N=2^14; % Aumentei para 11
t=0:h:(N-1)*h;

homg=2*pi/(N*h); % = Delta omega
% Com shift, os valores omega são:
omgmax=(N/2-1)*homg;
omgmin=(-N/2)*homg;
omg=omgmin:homg:omgmax;

nyquist=pi/h % Frequencia de Nyquist

%y = sin(10*t) + sin(10.05*t) 

y = sin(10*t) + sin(10.05*t);
z= fftshift(fft(y)); 
dens_esp=(h*abs(z)).^2;

figure(2)
subplot(2,1,1)
plot(omg, dens_esp)
title("\it{y}=sin(10t) + sin(10.05\it{t})")

subplot(2,1,2)
plot(t, y)

figure(3)
plot(omg, dens_esp)
axis([9.8 10.2 0 Inf])