%% Problema 5.1 - Condução de calor
% Feito com Euler
clc;
clear all;
close all;

alpha= 0.01;
u = 0.1;    

tf=50;
L=1; % [cm]

%Alinea b
dx=0.05;
dt=0.125;

% Alinea c
dx1=0.1;
dx2=0.01;

eta1= (2*alpha) / u; % critério de estabilidade
eta2= (dx^2)/(2*alpha);

% alinea c) Verificar se dx e dt satisfazem o critério de estabilidade
if (dx > eta1 && dt > eta2)
  error(['ERRO: Os valores de dt e dx não satisfazem o critério de ' ...
      'convergencia:\neta1= %.2f e eta2= %.2f'],eta1,eta2)
end

t=0:dt:tf;
x=0:dx:L;
nx=numel(x);
nt=numel(t);

T=zeros(nx,nt);
T(:,1)=100 *(x/L);
T(1,:)= 0; % 0 ºC
T(nx,:)=100; % 100 ºC

C = (u/2) * (dt/dx);
D = alpha * (dt/ dx^2); 

for n=1:nt-1
    for i=2:nx-1
        T(i,n+1)=T(i,n)+D*(T(i+1,n)-2*T(i,n)+T(i-1,n))- C*(T(i+1,n)-T(i-1,n));
    end
end

figure(1);
mesh(t(1:10:end),x,T(:,1:10:end))
xlim([0 tf])
xlabel('Tempo (s)')
ylim([0 L])
ylabel('Posicao x (cm)')
zlim([0 100])
zlabel('Temperatura ºC')

figure(2)

Test=100*(exp(10*x/L)-1)/(exp(10)-1);

plot(x,Test)
xlabel('x')
ylabel('T estacionário')

%%  CRITERIO DE ESTABILIDADE

nx1= 2*alpha/u;
nt1= dx1^2/(2*alpha);
disp([num2str(dx1), ' <= ', num2str(nx1),' e ' ,num2str(dt), ' <= ',num2str(nt1)])
disp('Logo dx=0.1 cm é estável')

nt2= dx2^2/(2*alpha);
disp([num2str(dx2), ' <= ', num2str(nx1),' e ' ,num2str(dt), ' <= ',num2str(nt2)])
disp('Logo dx=0.01 cm não é estável')