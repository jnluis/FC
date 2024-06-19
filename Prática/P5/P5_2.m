%% Problema 5.2 - Condução de calor
% Feito com Método de Crank-Nicolson
clc;
clear all;
close all;

kk= 0.93; % é o k dos dados do problema [cal/(s cm ºC)]
c=0.094; % [cal/(g ºC)]
ro=8.9;  % [g/cm^3]
D=kk/(c*ro); % fator k/(c*ro)

tf=500;
L=50; % [cm]
dx=0.5;
dt=0.1;

eta= D*dt/dx^2; % critério de estabilidade

t=0:dt:tf;
x=0:dx:L;
nx=numel(x);
nt=numel(t);

T=zeros(nx,nt);
T(1,:)= 0;
T(nx,:)=0;
T(2:nx-1,1)=100;

TT=T;
n_mat=nx-2;
aux_a=2/eta+2; % vem da solução
aux_b=2/eta-2;

alinea = input('Qual alinea?:\n a,b, ou c?\n','s');

% Matriz A
diagonal =aux_a*ones(n_mat,1);
A=diag(diagonal)-diag(ones(n_mat-1,1),-1)-diag(ones(n_mat-1,1),+1);

% [ aux_a   -1      0    ....   0 ]
% [   -1   aux_a   -1    ....   0 ]
% [   0     -1    aux_a  ....   0 ]
% [ ....    ....  ....   ....  ... ]
% [   0      0      0    ....  aux_a ]

%Vetor b
b= zeros(nx-2,1); % Para ser coluna

tic

for n=1:nt-1
    b=T(1:nx-2,n)+ aux_b*T(2:nx-1,n)+T(3:nx,n);
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
ylabel('Tempo (s)')

figure(2)
mesh(t,x,T)
xlim([0 tf])
xlabel('Tempo (s)')
ylim([0 L])
ylabel('Posicao x (cm)')
zlim([0 100])
zlabel('Temperatura ºC')