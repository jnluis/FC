%% Ex2a
clc;
close all;
clear all;

h0 = 0.35; 
g = 9.81;
D= 0.35;
d= 0.025;

t0= 0;
tf= 30;
h= 0.1;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = h0;

v = zeros(N,1);
v(1) = 0;

fx= @(V) V; % funções anónimas
fv = @(X) -sqrt(2*g)*((d/D).^2) * sqrt(X);

for k=1:N-1
    % Runge-Kutta           
    r1v = fv(x(k));
    r1x = fx(v(k));

    r2v = fv(x(k) +r1x*h/2);
    r2x = fx(v(k) +r1v*h/2);

    r3v = fv(x(k) +3*r2x*h/4);
    r3x = fx(v(k) +3*r2v*h/4);

    v(k+1) = v(k) + (2*r1v+3*r2v+4*r3v)*h/9;
    x(k+1) = x(k) + (2*r1x+3*r2x+4*r3x)*h/9;   
end

hA = (-sqrt(g/2)*((d/D).^2) *t + sqrt(h0)).^2;

figure(1)
subplot(1,2,1)
plot(t,x,t,hA)
title('RK3- Reservatório de água');
xlabel("Tempo (s)")
ylabel("Altura da superfície da água (m)")
