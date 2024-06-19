%*************************************************************************
%
% NOME 1: Bruno Figueiredo
% MEC 1: 103489
% Turma: PL8
% 
%*************************************************************************
%
% NOME 2: Laura Villalba
% MEC 2: 102847
% Turma: PL8
% 
%*************************************************************************
%
% NOME 2: Rafael Morgado
% MEC 2: 104277
% Turma: PL8
%
%*****************************Alínea b)***********************************
clc,clear,close all

%%  Variáveis

n=3;
w=3;
lameda=0.005;

h = 0.01; 
x = -3:h:3;
N = length(x);

%Potencial
V = 0.5*w^2*x.^2+lameda*x.^4;

y = zeros(1,N); 
dy = y; 
y(1) = 0; 
dy(1) = 0.01;


tol = 10^-10; 
result = zeros(1,100);
%Estimar os valores para o shooting
Enharm=(n+1/2)*w;   %Componente harmonica
%condiçoes para o shooting  (10.5 valor base para n=3)
E(1)=10.4;
E(2)=10.6;

%%  Método de shooting e Runge-Kutta de 4ª ordem
tic
for i = 1:1000
    % Funções anónimas

    fy = @(dy) dy; % "v"
    fDy = @(V,y) 2 * (V - E(i)) * y; % "f"

    for k = 1:N-1
        r1v = fDy(V(k),y(k));
        r1x = fy(dy(k));

        r2v = fDy(V(k),y(k) + r1x * h/2);
        r2x = fy(dy(k) + r1v * h/2);

        r3v = fDy(V(k),y(k) + r2x * h/2);
        r3x = fy(dy(k) + r2v * h/2);

        r4v = fDy(V(k),y(k) + r3x * h);
        r4x = fy(dy(k) + r3v * h);

        dy(k+1) = dy(k) + 1/6 * (r1v + 2 * r2v + 2 * r3v + r4v) * h;
        y(k+1) = y(k) + 1/6 * (r1x + 2 * r2x + 2 * r3x + r4x) * h;
    end

    % Inicio do "shoting"

    result(i) = y(end);

    if i > 1
        m = (result(i) - result(i-1))/(E(i) - E(i-1));
        E(i+1) = E(i) + (0 - result(i))/(m);
    end

    % Condição de paragem

    erro = abs(result(i)  - 0);

    if erro < tol
        break
    end
end
tempo=toc

%%  Gráfico e impressão do valor próprio

c = trapz(x,y.^2);
ynorm = y/sqrt(c);
plot(x,ynorm)
xlabel('x(m)')
ylabel('\psi_{norm}')
title("\psi_{norm} em função de x")
grid on

fprintf('E3 = %f \n',E(end)) %Valor próprio da energia 