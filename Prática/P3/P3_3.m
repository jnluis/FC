%% Problema 3.3 - Sistema Massa/Mola RK de passo adaptativo
clc;
close all;
clear all;

x0 = 1; % posição inicial
v0 = 0;

K  = 16; %N/m
m = 1; % 1 kg

w = sqrt(K/m);
w2 = K/m;

t0= 0;
tf= 10;
h=0.01;
t = t0:h:tf;

options=odeset('RelTol',3e-14, "AbsTol",[1e-13 1e-13]);
[t, solucao] = ode45(@f3_3, 0:h:tf, [x0 v0], options, m, K);

% RelTol: Define a tolerância relativa para a solução
% Valores menores de RelTol resultarão numa solução mais precisa, 
% mas também pode aumentar o tempo de execução

% AbsTol: Define a tolerância absoluta para a solução
% È um vetor contendo uma tolerância absoluta para cada c
% [abstol_1 abstol_2] indica tolerâncias absolutas separa

x = solucao(:,1);
v = solucao(:,2);

n= numel(t);

E = 0.5*K*x.^2+0.5*m*v.^2; % Runge-Kutta

figure(1)
plot(t,x)
xlabel('t');ylabel('x');

figure(2)
plot(x,v)
xlabel('x');ylabel('v');

figure(3)
plot(t,E)
xlabel('t');ylabel('E');