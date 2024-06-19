% Parte 2
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
%*****************************Alínea b)***********************************
clc;
clear all;
close all;

a = [0 2 3 5 6]; % Diferentes valores de alfa

N = 1024;
L = 80;
dx = L/N;
x = -(N-1)/2*dx : dx : (N-1)/2*dx;

ti = 0;
tf = 1;
dt = 1.5e-5;
t = ti : dt : tf;
Nt = length(t);

% Vetor das frequências angulares
dw = 2*pi / (N*dx);
wmax = (N/2-1);
wmin = (-N/2);
w = [0:wmax wmin:-1] * dw;
 
qx = zeros(Nt, N);
qx(1, :) = -12.*sech(x).^2;
t1 = (1i.*w);
t3 = (1i.*w).^3;

alphas = [0, 2, 3, 5, 6];
colors = ['b','g','r','c','m'];
figure;

for is = 1:length(a)
    alfa = a(is);
    
    for n = 1:Nt-1
        q = qx(n, :);
        r1 = (-ifft(t3.*fft(q)) + alfa.*q.*ifft(t1.*fft(q)));
        v = q + r1*dt/2;
        r2 = (-ifft(t3.*fft(v)) + alfa.*v.*ifft(t1.*fft(v)));
        v2 = q + r2*dt/2;
        r3 = (-ifft(t3.*fft(v2)) + alfa.*v2.*ifft(t1.*fft(v2)));
        v3 = q + r3*dt;
        r4 = (-ifft(t3.*fft(v3)) + alfa.*v3.*ifft(t1.*fft(v3)));

        qx(n+1, :) = qx(n, :) + (r1 + 2*r2 + 2*r3 + r4)*dt/6;
    end
    intensidade = qx(end, :).^2;
    plot(x, intensidade(end,:), colors(is), 'DisplayName', ['\alpha = ' num2str(alfa)]);
    hold on
end
title('Perfis finais em função de x para diferentes valores de \alpha');
xlabel('x');
ylabel('Intensidade');
legend;
grid on;
