clc;
clear all;
close all;

g = 9.8;
v0 = 0;
vlim = 6.8;
tf = 0.5; % Tempo final modificado para 0.5 s

h_values = 10.^(-3:-1:-7); % Valores de h
errors = zeros(size(h_values)); % Para armazenar os erros para cada h

for i = 1:length(h_values)
    h = h_values(i);
    t = 0:h:tf;
    N = length(t);
    v = zeros(N,1);
    v(1) = v0; % Velocidade inicial
    
    % Cálculo numérico da velocidade
    for k = 1:N-1
        v(k+1) = v(k) + (((-g * abs(v(k))* v(k)) / vlim^2) - g) *h;
    end
    
    % Valor analítico de vz em t=0.5 s
    vanali = -vlim * tanh((g/vlim)*tf);
    
    % Cálculo do erro
    errors(i) = abs(vanali - v(end));
end

% Plot do logaritmo do erro em função do logaritmo de h
log_h = log(h_values);
log_errors = log(errors);

figure(1);
plot(log_h, log_errors, 'o-');
title('Logaritmo do Erro vs. Logaritmo de h');
xlabel('Logaritmo de h');
ylabel('Logaritmo do Erro');
grid on;
lsline; % Linha de reta média

% Determinação do expoente usando polyfit
coefficients = polyfit(log_h, log_errors, 1);
fprintf("O expoente aproximado é %d \n", coefficients(1));
