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

y = sin(t);
z= fftshift(fft(y));
dens_esp=(h*abs(z)).^2;

% 𝑦 = sin(𝑡)
subplot(3,2,1)
plot(omg, dens_esp)
title("\it{y}=sin(\it{t})")

subplot(3,2,2)
plot(t, y)

%y = 10 sin(𝑡)

y = sin(10*t);
z= fftshift(fft(y));
dens_esp=(h*abs(z)).^2;

subplot(3,2,3)
plot(omg, dens_esp)
title("\it{y}=sin(10\it{t})")

subplot(3,2,4)
plot(t, y)

%y = sin(t) + 10 sin(𝑡)

y = sin(t) + sin(10*t);
z= fftshift(fft(y)); 
dens_esp=(h*abs(z)).^2;

subplot(3,2,5)
plot(omg, dens_esp)
title("\it{y}=sin(t) + sin(10\it{t})")

subplot(3,2,6)
plot(t, y)