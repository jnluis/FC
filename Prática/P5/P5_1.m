%% Problema 5.1 - Condução de calor
% Feito com Euler
clc;
clear all;
close all;

kk= 0.93; % é o k dos dados do problema [cal/(s cm ºC)]
c=0.094; % [cal/(g ºC)]
ro=8.9;  % [g/cm^3]
D=kk/(c*ro);  % fator k/(c*ro)

tf=500;
L=50; % [cm]
dx=0.5;
dt=0.1;

eta= D*dt/dx^2; % critério de estabilidade

% 5.1b) Verificar se dx e dt satisfazem o critério de estabilidade
if (eta > 0.5)
  error(['ERRO: Os valores de dt e dx não satisfazem o critério de ' ...
      'convergencia:\neta= %.2f'],eta)
end

t=0:dt:tf;
x=0:dx:L;
nx=numel(x);
nt=numel(t);

T=zeros(nx,nt);
T(1,:)= 0;
T(nx,:)=0;
T(2:nx-1,1)=100;

for n=1:nt-1
    for i=2:nx-1
        T(i,n+1)=T(i,n)+eta*(T(i+1,n)-2*T(i,n)+T(i-1,n));
    end
end

figure(1);
contourf(x,t,T')
colorbar
xlabel('Posicao x (cm)')
ylabel('Tempo (s)')

figure(2)
mesh(t,x,T)
xlim([0 tf])
xlabel('Tempo (s)')
ylim([0 L])
ylabel('Posicao x (cm)')
zlim([0 100])
zlabel('Temperatura ºC')

%% Alínea C

i_quarto=(nx-1)/4+1;
T_quarto=interp1(T(i_quarto,i_quarto:end),t(i_quarto:end),T(i_quarto,2)/2,'linear');
fprintf('Tempo: %f \n', T_quarto)
figure(3)
plot(t,T(i_quarto,:));