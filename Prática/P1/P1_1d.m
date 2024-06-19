%% Problema 1.1 - Queda de uma Pedra
clc;
clear all;
close all;

t0 =0;
tf= 2;
h = 0.2; % passo
t = t0:h:tf;

g= 9.8;
m = 0.150; % já em kg
v0=0;
H= 6;

N = length(t);
v = zeros(N,1);
z = zeros(N,1);
v(1) = v0; % a velocidade inicial para este problema é 0
z(1) = H; % porque é a altura de 2 andares

for k=1:N-1
    v(k+1) = v(k) -g *h;
    z(k+1) = z(k) + v(k)*h;

end

vanali = v0 - g *t;
zanali = H - 0.5*g * (t.^2);

%subplot(1, 2, 1);
figure;
plot(t,z,"b");
hold on;
plot(t,zanali, 'r');
legend("numérico","analítico")
hold off;
xlim([0 1.2]);
title('Posições ao longo do tempo');
xlabel("Tempo (s)")
ylabel("Posição (m)")

%subplot(1, 2, 2);
figure;
plot(t,v);
hold on
plot(t, vanali);
xlim([0 1.2]);
title('Velocidade ao longo do tempo')
legend("Numérico","Analítico");
xlabel("Tempo (s)")
ylabel("Velocidade (m/s)")

% Encontre o primeiro índice onde z <= 0
idx_solo = find(z <= 0, 1);

% Instante em que cai no chão
t_solo = interp1([z(idx_solo-1), z(idx_solo)], [t(idx_solo-1), t(idx_solo)], 0, "linear");
v_solo = interp1([z(idx_solo-1), z(idx_solo)], [v(idx_solo-1), v(idx_solo)], 0, "linear");


fprintf("Instante em que chega ao solo = %f s\n\n", t_solo)
fprintf("Velocidade com que chega ao solo = %f m/s\n\n", v_solo)