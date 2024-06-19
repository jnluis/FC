% Parte 1
%*************************************************************************
%
% NOME 1: Martim Gil
% MEC  1: 102901
% Turma : PL6
% 
%*************************************************************************
%
% NOME 2: João Luís
% MEC  2: 107403
% Turma : PL6
% 
%*************************************************************************
%
% NOME 3: João Marques
% MEC  3: 108072
% Turma : PL6
%
%*************************************************************************
%
% NOME 4: Délcio Amorim
% MEC  4: 109680
% Turma : PL6
%
%*****************************Alínea a)***********************************
clc;
clear all;
close all;

% Parâmetros iniciais
h = 0.05;
N = 1024;

% Variar o h e o N para ver o que acontece
%h=0.3;
%N=512;
x = -N/2*h:h:N/2*h; %domínio tem de ser simétrico em 0

q = tanh(x).*sech(x);
% Transformada de Fourier
Q = fftshift(fft(ifftshift(q)));

% Frequências associadas aos pontos no domínio da frequência
f = (-N/2:N/2)/(N*h);

% Derivação no domínio da frequência
D4Q = (1i*2*pi*f).^4 .* Q;
% Transformada de Fourier inversa para obter a derivada
d4q = fftshift(ifft(ifftshift(D4Q)));

% Quarta derivada analítica que vem do Wolfram Alpha
d4q_analytical = tanh(x) .* sech(x) .* (tanh(x).^4 + 61*sech(x).^4 - 58*tanh(x).^2 .* sech(x).^2); 

figure(1);
grid on
plot(x, d4q, 'r--', x, d4q_analytical, 'k');
legend('Derivada numérica', 'Derivada analítica');
title('Comparação entre a quarta derivada numérica e analítica');
xlim([-10 10])
xlabel('x');
ylabel('d4q/dx^4');