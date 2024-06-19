%exercicio 3)
clc;
clear all;
close all;

%constantes
%h pode variar
h=0.02;
t=0:h:20;
t0=0;
tf=50;
x0=1;
v0=1;

N=length(t);

alfa=-0.1;
m=1;
K=1;

x=zeros(1,N);
x(1)=1;

vx = zeros(1,N);
vx(1)= 1;

reltol = 3e-14;
abstol_1 = 1e-13;
abstol_2 = 1e-13;

options = odeset('RelTol', reltol, 'AbsTol', [abstol_1 abstol_2]);
[t, sol] = ode45(@(t, sol)f(t, sol, K, m), [t0 tf], [x0 v0], options);

x = sol(:, 1);
v = sol(:, 2);

plot(t, x);
title('Posição em função do tempo de um sistema massa-mola, usando a função ode45');
xlabel('t / s');
ylabel('x / m');