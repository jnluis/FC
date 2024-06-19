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
vx0=0;

N = numel(t);
x = zeros(N,1);
x(1) = h0;

v = zeros(N,1);
v(1) = vx0;

fx= @(V) V; % funções anónimas
fv = @(X) -sqrt(2*g)*((d/D).^2) * sqrt(X);

m=1;
K=1;
options=odeset('RelTol',3e-14, "AbsTol",[1e-13 1e-13]);
[t, solucao] = ode45(@f, t0:h:tf, [h0 vx0], options, m, K);

x = solucao(:,1);
v = solucao(:,2);

n= numel(t);

hA = (-sqrt(g/2)*((d/D).^2) *t + sqrt(h0)).^2;

figure(1)
subplot(1,2,1)
plot(t,x,t,hA)
title('RK3- Reservatório de água');
xlabel("Tempo (s)")
ylabel("Altura da superfície da água (m)")
