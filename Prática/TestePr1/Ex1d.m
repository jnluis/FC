%% Ex1d
clc;
close all;
clear all;

x0 = 0.1; % em m
vx0 = 0;
y0= 0.1;
vy0 = 0;

K  = 0.4; %N/m
m = 0.2; % 200 g
L = 1;

w = K/m;
w2 = w^2;
g = 9.81;

t0= 0;
tf= 30;
h= 0.0001;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

vx = zeros(N,1);
vx(1) = vx0;

y = zeros(N,1);
y(1) = y0;

vy = zeros(N,1);
vy(1) = vy0;

for k=1:N-1
    vx(k+1) = vx(k) + (-g*(x(k)/L) - K*(x(k)-y(k)))*h;
    vy(k+1) = vy(k) + (-g*(y(k)/L) + K*(x(k)-y(k)))*h;
    
    x(k+1) = x(k) + vx(k+1) *h;
    y(k+1) = y(k) + vy(k+1) *h;
end

w1= sqrt(g/L);
w2 = sqrt((g/L)+2*(K/m));
xa= (0.5 * (x0 +y0)* cos(w1*t) ) + (0.5 * (x0 - y0)* cos(w2*t) ); 

figure(1)
plot(t,x)
hold on
plot(t,xa)
title('X ao longo do tempo');
legend('numerico','analitico');
xlabel("Tempo (s)")
ylabel("X (m)")
hold off

figure(2)
plot(t,y,t,xa)
legend('numerico','analitico');
title('Y ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Y (m)")


% plot(t,vx)
% title('Velocidades ao longo do tempo');
% xlabel("Tempo (s)")
% ylabel("Velocidade (m/s)")
% 
% subplot(2,2,3)
% plot(x,vx)
% xlabel("Posição (m)")
% ylabel("Velocidade (m/s)")
% 
% subplot(2,2,4)
% plot(t,Em)
% title('Energia mecânica ao longo do tempo');
% xlabel("Tempo (s)")
% ylabel("Em (J)")

