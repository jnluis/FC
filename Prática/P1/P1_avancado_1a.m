%% Problema 1.1 - Movimento a uma dimensão de um volante de badminton
clc;
clear all;
close all;

t0 =0;
tf= 0.5;
h = 0.1; % passo
t = t0:h:tf;

g= 9.8;
v0=0;
vlim= 6.8; 

N = length(t);
v = zeros(N,1);
z = zeros(N,1);
a = zeros(N,1);
v(1) = v0; % a velocidade inicial para este problema é 0

for k=1:N-1
    v(k+1) = v(k) + (((-g * abs(v(k))* v(k)) / vlim^2) - g) *h;
end

vanali = -vlim * tanh((g/vlim)*t);

%subplot(1, 2, 2);
figure;
plot(t,v);
hold on
plot(t, vanali);
title('Velocidade ao longo do tempo')
legend("Numérico","Analítico");
xlabel("Tempo (s)")
ylabel("Velocidade (m/s)")