% Parte 1
%*************************************************************************
%
% NOME 1: Martim Gil
% MEC  1: 102901
% Turma : PL6
% 
%*************************************************************************
%
% NOME 2: João Luís
% MEC  2: 107403
% Turma : PL6
% 
%*************************************************************************
%
% NOME 3: João Marques
% MEC  3: 108072
% Turma : PL6
%
%*************************************************************************
%
% NOME 4: Délcio Amorim
% MEC  4: 109680
% Turma : PL6
%
%*****************************Alínea b)***********************************
clc;
clear all;
close all;

% Parâmetros iniciais
h = 0.05;
N = 1024;
x = -N/2*h:h:N/2*h; %domínio tem de ser simétrico em 0

q = tanh(x).*sech(x);

% Para diferenças finitas, precisamos do número de pontos
N = length(x);
    
% Inicializa o vetor de derivadas
d4q = zeros(1, N);
    
for k = 3:(N-2)
    d4q(k) = (q(k-2) - 4*q(k-1) + 6*q(k) - 4*q(k+1) + q(k+2)) / h^4;
end
    
% Ajusta os pontos nas bordas (não calculados pelo esquema central)
% Pode-se usar esquemas de ordem menor ou extrapolação, aqui vamos deixar como zero
%y_deriv(1:2) = NaN; % Ou outro valor conforme a necessidade
%y_deriv(N-1:N) = NaN; % Ou outro valor conforme a necessidade

% Quarta derivada analítica que vem do Wolfram Alpha
d4q_analytical = tanh(x) .* sech(x) .* (tanh(x).^4 + 61*sech(x).^4 - 58*tanh(x).^2 .* sech(x).^2); 

figure(1);
grid on
plot(x, d4q, 'r--', x, d4q_analytical, 'k');
legend('Derivada numérica', 'Derivada analítica');
title('Comparação entre a quarta derivada numérica e analítica');
xlim([-10 10])
xlabel('x');
ylabel('d4q/dx^4');