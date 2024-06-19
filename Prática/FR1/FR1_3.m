%% Problema 2.3 - Oscilador Quártico - ode45
clc;
close all;
clear all;

x0 = 1; % posição inicial
vx0 = 1;
alpha = -0.1;

K  = 1; %N/m
m = 1; % 1 kg

w = sqrt(K/m);
w2 = K/m;

t0= 0;
tf= 50;
h= 0.01;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

v = zeros(N,1);
v(1) = vx0;

% Função para a EDO
%fx = @(t,x,vx) vx;
%fv = @(t,x,vx) -1 * (x+2 * alpha * x^3);

%[t, y] = ode45(@(t,y) [fx(t,y(1), y(2)); fv(t, y(1), y(2))], t0:h:tf, [x0 vx0]);
% x = y(:,1);
% vx = y(:,2);

options=odeset('RelTol',3e-14, "AbsTol",[1e-13 1e-13]);
[t, solucao] = ode45(@f, t0:h:tf, [x0 vx0], options, m, K);

% RelTol: Define a tolerância relativa para a solução
% Valores menores de RelTol resultarão numa solução mais precisa, 
% mas também pode aumentar o tempo de execução

% AbsTol: Define a tolerância absoluta para a solução
% È um vetor contendo uma tolerância absoluta para cada c
% [abstol_1 abstol_2] indica tolerâncias absolutas separa

x = solucao(:,1);
v = solucao(:,2);

n= numel(t);

Em = 0.5*(K*(x.^2+alpha*x.^4)+ m*v.^2);

figure(2)
plot(x,v)
title('ode45- Oscilador Quártico');
xlabel('x');ylabel('v');

figure(3)
plot(t,Em)
title('ode45- Oscilador Quártico');
xlabel('t');ylabel('Em');