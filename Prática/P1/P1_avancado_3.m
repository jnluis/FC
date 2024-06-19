%% Problema 1.3 Bola de futebol — desvios laterais
clc;
close all;
clear all;

t0= 0;
tf= 2;
h= 0.001;
t = t0:h:tf;
N = length(t);
g = 9.8;

m = 0.45;
v0 = 80 /3.6; % Para ficar em m/s
alpha = 10; % em graus
H = 0; % Altura

r=600; % em rpm, precisamos de w em rad/s
w = (pi *r) / 30; % velocidade angular

p= 0.7; % é o perímetro, precisamos da área
raio=p/(2*pi); %raio da bola

A = pi * (raio^2);
Cm =1;

c= 0.5*Cm * A*1.225*raio; % é o valor da Força de Magnus excetuando w * v

r= zeros(N,3); % r é x, y e z
v_vec = zeros(N,3);

a = [0 0 -g];
r(1,:) = [0 0 H];
v_vec(1,:)=[v0*cosd(alpha) 0 v0*sind(alpha)]; % o d é porque o ângulo está em graus
w_vec = [0 0 w]; % só no eixo do Z?

for k=1:N-1
    v=norm(v_vec(k,:));
    if v<=9
        Fd=-(0.015*v^2)*v_vec(k,:)/v;
    elseif v>9 & v<=20
        Fd=-(0.25147+0.17431*v-0.015*v^2+0.00054*v^3)*v_vec(k,:)/v; % O 0.015 está a fazer com que dê diferente das soluções!
    % ela tem este valor, 0.01384, que é do antigo
    else
        Fd=-(-4.025+0.323*v)*v_vec(k,:)/v;
    end
    Fl=c*cross(w_vec,v_vec(k,:)); % é a força de magnus. O cross faz o produto externo
    a=(Fd+Fl)/m+[0 0 -g];

    v_vec(k+1,:)=v_vec(k,:)+a*h;
    r(k+1,:)=r(k,:)+v_vec(k,:)*h; 

    if r(k+1,3)<0 % quando a bola volta ao chão 
        break
    end
end

t=t(1:k+1);
r=r(1:k+1,:);
v_vec=v_vec(1:k+1,:);

figure;
plot3(r(:,1),r(:,2), r(:,3));
%xlim([0 1.2]);
title('Trajetória da bola com desvio lateral')
legend("trajetória");
xlabel("x")
ylabel("y")
zlabel("z")

%% Podemos interpolar para encontrar o instante em que bateu no chão 
%% e a partir daí calcular a distância percorrida

desvio=interp1(r(end-1:end,3),r(end-1:end,2),0,'linear')
% Calcula a distância horizontal até o ponto de impacto no solo
distancia_horizontal_exata = interp1(r(k:k+1,3), r(k:k+1,1), 0, 'linear')