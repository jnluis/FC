%% Problema 5.3 - Condução de calor
% Feito com Método de Crank-Nicolson
clc;
clear all;
close all;

kk= 0.93; % é o k dos dados do problema [cal/(s cm ºC)]
c1=0.094; % [cal/(g ºC)]
c2=0.188;

ro=8.9;  % [g/cm^3]
D1=kk/(c1*ro); % fator k/(c*ro)
D2=kk/(c2*ro); % fator k/(c*ro)

tf=500;
L=50; % [cm]
dx=0.5;
dt=0.1;

eta1= D1*dt/dx^2; % critério de estabilidade
eta2= D2*dt/dx^2; % critério de estabilidade

t=0:dt:tf;
x=0:dx:L;
nx=numel(x);
nt=numel(t);

T=zeros(nx,nt);
T(1,:)= 0;
T(nx,:)=0;
T(2:nx-1,1)=50*sin(2*pi*x(2:nx-1)/L); % Temperatura inicial

TT=T;
n_mat=nx-2;
n_mat_meio=(n_mat+1)/2; % ponto central

% Matriz A
diagonal(1:n_mat_meio-1)=2/eta1+2;
diagonal(n_mat_meio)=4/(eta1+eta2)+2;
diagonal(n_mat_meio+1:n_mat)=2/eta2+2;
A=diag(diagonal)-diag(ones(n_mat-1,1),-1)-diag(ones(n_mat-1,1),+1);

% 2/eta_variavel-2 para escrita do b
aux_b=zeros(nx-2,1);
aux_b(1:n_mat_meio-1)=2/eta1-2;
aux_b(n_mat_meio)=4/(eta1+eta2)-2;
aux_b(n_mat_meio+1:n_mat)=2/eta2-2;

b=zeros(nx-2,1); % Para ser coluna

alinea = input('Qual alinea?:\n a,b, ou c?\n','s');

tic
for n=1:nt-1
    b=T(1:nx-2,n)+ aux_b.*T(2:nx-1,n)+T(3:nx,n);
    b(1)= b(1)+T(1,n+1);
    b(n_mat)=b(n_mat)+T(nx,n+1);

    switch alinea
        case 'a'
            T(2:nx-1,n+1)= linsolve(A,b);
        case 'b'
            T(2:nx-1,n+1)=sol_sist_trid(A,b);
        case 'c'
            [LL,U,P] = lu(A);
            y=LL\b;
            T(2:nx-1,n+1)=U\y;
        otherwise
            fprintf('Nope\n\n')
            return
    end

end
toc

figure(1);
contourf(x,t,T')
colorbar
xlabel('Posicao x (cm)')
ylim([0 250])
ylabel('Tempo (s)')

figure(2)
mesh(t,x,T)
xlim([0 tf])
xlabel('Tempo (s)')
ylim([0 L])
ylabel('Posicao x (cm)')
zlim([-50 50])
zlabel('Temperatura ºC')