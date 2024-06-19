%% Problema 3.2 - Sistema Massa/Mola com RK4
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

Hhs = 10.^(-1:-1:-6); % Passo h
NH = length(Hhs); % Nº de pontos associados ao passo h

ERRO = zeros(1,NH); % Inicializa o vetor de erros

for j = 1:NH
    h = Hhs(j);
    t = t0:h:tf;
    N = numel(t);

    x = zeros(N,1);
    x(1) = x0;

    v = zeros(N,1);
    v(1) = v0;

    fx= @(V) V; % função anónima
    fv = @(X) -K*X/m;
    for k = 1:N-1
        r1v = fv(x(k));
        r1x = fx(v(k));
        r2v = fv(x(k) +r1x*h/2);
        r2x = fx(v(k) +r1v*h/2);
        r3v = fv(x(k) +r2x*h/2);
        r3x = fx(v(k) +r2v*h/2);
        r4v = fv(x(k) +r3x*h);
        r4x = fx(v(k) +r3v*h);

        v(k+1) = v(k) + (r1v+2*r2v+2*r3v+r4v)*h/6;
        x(k+1) = x(k) + (r1x+2*r2x+2*r3x+r4x)*h/6;
    end

    % Cálculo do erro no último ponto
    ERRO(j) = abs(x(j+1) - x(j));
    logERRO = log(ERRO);

end

logHhs = log(Hhs);

figure(1);
plot(logHhs, logERRO, '+');
xlabel("ln (h)")
ylabel("ln (Erro)")
legend("Euler", "Location", "north")
lsline

% Determinação do expoente usando polyfit
coefficientsE = polyfit(logHhs, logERRO, 1);
fprintf("Slope RK4 %d \n", coefficientsE(1));