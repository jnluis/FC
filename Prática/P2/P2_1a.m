%% Problema 2.1 - Sistema Massa/Mola
clc;
close all;
clear all;

x0 = 1; % posição inicial
vx0 = 0;

K  = 1; %N/m
m = 1; % 1 kg

w = K/m;
w2 = w^2;

t0= 0;
tf= 50;
h= 0.01;
t = t0:h:tf;

N = numel(t);
x = zeros(N,1);
x(1) = x0;

vx = zeros(N,1);
vx(1) = vx0;

n = input(['Escolha o método: \n' ...
    '1- Euler; \n' ...
    '2- Euler-Cromer;\n' ...
    '3- Euler implícito sem linsolve;\n' ...
    '4- Euler implícito com linsolve;\n' ...
    '5- Crank-Nicolson com linsolve;\n']);

switch n
    case 1
        disp("Método de Euler")
        for k=1:N-1
            a = -K* x(k) /m;
            vx(k+1) = vx(k) +a *h;
            x(k+1) = x(k) + vx(k) *h;
        end
    case 2
        disp("Método de Euler-Cromer")
        for k=1:N-1
            a = -K* x(k) /m;
            vx(k+1) = vx(k) +a *h;
            x(k+1) = x(k) + vx(k+1) *h;
        end
    case 3
        disp("Método de Euler Implícito sem linsolve")
        aux= 1 + w *h^2; % w2= K/m
        for k=1:N-1
            
            % Atualização da posição usando o método de Euler implícito
            x(k+1) = (x(k) + vx(k) *h) / aux;

            % Atualização da velocidade
            vx(k+1) = vx(k) -w2*x(k+1)*h;
            
        end
    case 4
        disp("Método de Euler Implícito com linsolve")
        A = [1 -h; w^2*h 1];
        for k=1: N-1
            % Definindo o vetor b com base nas equasções do método de Euler
            % Implícito
            b = [x(k); vx(k)];

            % Resolvendo o sistema linear utilizando linsolve
            % X = linsolve( A, B) solves the linear system AX = B
            aux = linsolve(A, b);

            % Atualizando as variáveis de posição e velocidade
            vx(k+1) = aux(2);
            x(k+1) = aux(1);
        end
    case 5
        disp("Método de Crank-Nicolson com linsolve")
        A = [1 -h/2; w^2*h/2 1];
            for k=1: N-1
            % Definindo o vetor b com base nas equasções do método de Euler
            % Implícito
            b = [x(k)+vx(k)*h/2; vx(k)-w^2*x(k)*h/2];
                 % posição         % velocidade

            % Resolvendo o sistema linear utilizando linsolve
            % X = linsolve( A, B) solves the linear system AX = B
            aux = linsolve(A, b);

            % Atualizando as variáveis de posição e velocidade
            vx(k+1) = aux(2);
            x(k+1) = aux(1);
            end
    otherwise
        fprintf('Nope');
        return
end
Ec = m*vx.^2;
Ep = K*x.^2;
Em = 0.5* (Ec + Ep);
subplot(2,2,1)
plot(t,x)
title('Posições ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Posição (m)")

subplot(2,2,2)
plot(t,vx)
title('Velocidades ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Velocidade (m/s)")

subplot(2,2,3)
plot(x,vx)
xlabel("Posição (m)")
ylabel("Velocidade (m/s)")

subplot(2,2,4)
plot(t,Em)
title('Energia mecânica ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Velocidade (m/s)")