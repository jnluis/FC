clc 
clear all
close all

%input que pode ser modificado
x0 = 1.0;
vx0 = 0.0;

%parametros dos metodos numericos 
h = input('Introduza o valor de h em segundos: ');
tfin = 100;

%constantes do problema 
K = 1; % K MAIUSCULO
m = 1;

%calculos auxiliares
w = sqrt(K/m); 
N = numel(t);
x = zeros(N,1);
x(1) = x0;
vx = zeros(N,1);
vx(1)=vx0;

%escolha do método

met = ["Euler", "Euler-Cromer", "Euler implícito sem linsolve", ...
    "Euler implícito com linsolve", "Crank-Nicolson com linsolve"];
promtn = sprintf(["De entre os métodos, 1, 2, 3, 4, 5, escolha um: "], met);
n= input(promptn);

switch n
    case 1
        for k = 1:N-1
            a = -K*x(k)/m;
            vx(k+1) = vx(k) * a*h;
            x(k+1) = x(k)+vx(k)*h;
        end

    case 2
        for k = 1:N-1
            a = -K*x(k)/m;
            vx(k+1) = vx(k) * a*h;
            x(k+1) = x(k)+vx(k)*h;
         end