clc;
close all;
clear all;

L = 0.25; % indutância
C = 1e-3; % capacidade
V0 = 5;

t0 = 0;
tf = 0.5;

a = 1/(L*C);
Hhs = 10.^(-3:-1:-7); % Passo h
NH = length(Hhs); % Nº de pontos associados ao passo h

ERRO = zeros(1,NH); % Inicializa o vetor de erros
ERRO_EC = zeros(1,NH);

for j = 1:NH
    h = Hhs(j);
    t = t0:h:tf;
    N = numel(t);

    W = sqrt(1/(L*C));

    % Solução analítica
    Vt = V0 * cos(W*t);
   
    VE = zeros(1,N); % VC
    dVE = zeros(1,N);
    VE(1) = V0;
    dVE(1) = 0;

    VEc = zeros(1,N);
    dVEc = zeros(1,N);
    VEc(1) = V0;
    dVEc(1) = 0;

    for k = 1:N-1
        dVE(k+1) = dVE(k) + (-a*VE(k)) * h; % EULER  
        dVEc(k+1) = dVEc(k) + (-a*VEc(k)) * h; % EULER-CROMER

        VE(k+1) = VE(k) + dVE(k) * h; % EULER  
        VEc(k+1) = VEc(k) + dVEc(k+1) * h; % EULER-CROMER
    end

    % Cálculo do erro no último ponto
    ERRO(j) = abs(VE(end) - Vt(end));
    logERRO = log(ERRO);

    ERRO_EC(j) = abs(VEc(end) - Vt(end));
    logERRO_EC = log(ERRO_EC);
end

logHhs = log(Hhs);

figure(1);
plot(logHhs, logERRO, 'o');
xlabel("log (h)")
ylabel("log (| Erro |)")
legend("Euler", "Location", "north")
lsline

figure(2);
plot(logHhs, logERRO_EC, 'o');
xlabel("log (h)")
ylabel("log (| Erro |)")
legend("Euler-Cromer", "Location", "north")
lsline

% Determinação do expoente usando polyfit
coefficientsE = polyfit(logHhs, logERRO, 1);
coefficientsEc = polyfit(logHhs, logERRO_EC, 1);
fprintf("Slope Euler %d \n", coefficientsE(1));
fprintf("Slope Euler-Cromer %d \n", coefficientsEc(1));