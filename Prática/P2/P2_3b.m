%% Problema 2.3 - Oscilador Quártico
%% Alínea b' - Como varia o período com a amplitude? Faça a amplitude variar de 0.1 a 2 m
%% e represente graficamente os resultados.
clc;
close all;
clear all;

% Parâmetros fixos
alpha = -0.1;
K  = 1; %N/m
m = 1; % 1 kg
w = K/m;
t0= 0;
tf= 50;
h= 0.01;

% Amplitudes iniciais
amplitudes = 0.1:0.1:2;
periodos = zeros(size(amplitudes));

for j = 1:length(amplitudes)
    x0 = amplitudes(j);
    vx0 = 0; % Velocidade inicial para máxima amplitude
    
    t = t0:h:tf;
    N = numel(t);
    x = zeros(N,1);
    x(1) = x0;
    vx = zeros(N,1);
    vx(1) = vx0;
    
    for k=1:N-1
        a = -K*( x(k) + 2* alpha* x(k).^3 )/m;
        vx(k+1) = vx(k) +a *h;
        x(k+1) = x(k) + vx(k+1) *h;
    end
    
    % Localiza máximos
    imax = 0;
    tmax = [];
    xmax = [];
    for k= 2:N-1
        if and(x(k+1)-x(k) <=0 , x(k) - x(k-1) >= 0 )
            imax = imax+1;
            aux = lagr(t(k-1:k+1), x(k-1:k+1));
            tmax(imax) = aux(1);
            xmax(imax) = aux(2);
        end
    end
    
    if imax > 1 % Garante que pelo menos um período foi encontrado
        plf = polyfit(1:imax, tmax, 1);
        T = plf(1);
    else
        T = NaN; % Não foi possível determinar o período
    end
    
    periodos(j) = T;
end

figure(2);
plot(amplitudes, periodos, '*');
xlabel('Amplitude Inicial (m)');
ylabel('Período (s)');
title('Variação do Período com a Amplitude');
xlim([0 2.5])