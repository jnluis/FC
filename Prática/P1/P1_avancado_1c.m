%% Problema 1.1 - Movimento a uma dimensão de um volante de badminton
%% Lançamento para cima
clc;
clear all;
close all;

t0 =0;
tf= 2.5;
h = 0.001; % passo
t = t0:h:tf;

g= 9.8;
v0=16; % lançado para cima
y0 = 1;
vlim= 6.8; 

N = length(t);
v = zeros(N,1);
z = zeros(N,1);
a = zeros(N,1);
z = zeros(N,1);
v(1) = v0; % a velocidade inicial para este problema é 16 m/s
z(1) = y0;

for k=1:N-1
    v(k+1) = v(k) + (((-g * abs(v(k))* v(k)) / vlim^2) - g) *h;
    z(k+1) = z(k) + v(k) * h; % posiçao no instante
    if z(k+1) <0
        break
    end
end
%% Deveria ser feito este reassign para só ter um valor negativo??
% z = z(1:k+1); % se colocar apenas k em vez de k+1, só vou ter valores
% positivos, e depois não dá para interpolar
% t = t(1:k+1);
% v = v(1:k+1);
%% Da maneira que está, chega ao instante negativo e depois é sempre 0
vanali = -vlim * tanh((g/vlim)*t); % Talvez este analítico não dê para aplicar aqui

%subplot(1, 2, 2);
figure(2);
plot(t,v);
hold on
plot(t, vanali);
title('Velocidade ao longo do tempo')
legend("Numérico","Analítico");
xlabel("Tempo (s)")
ylabel("Velocidade (m/s)")
hold off

figure(1);
plot(t,z);
title('Posição ao longo do tempo')
legend("Numérico");
xlabel("Tempo (s)")
ylabel("Posição (m)")

% Encontrar o primeiro índice onde z <= 0
idx_solo = find(z <= 0, 1);

% Instante em que cai no chão
t_solo = interp1([z(idx_solo-1), z(idx_solo)], [t(idx_solo-1), t(idx_solo)], 0, "linear");
v_solo = interp1([z(idx_solo-1), z(idx_solo)], [v(idx_solo-1), v(idx_solo)], 0, "linear");

fprintf("Instante em que chega ao solo = %f s\n\n", t_solo)
fprintf("Velocidade com que chega ao solo = %f m/s\n\n", v_solo)