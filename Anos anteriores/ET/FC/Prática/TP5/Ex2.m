%%  Ex5.2

close,clc,clear all

%% Constantes

L=50;   %cm

k=0.93;
c=0.094;
p=8.9;

tf=500;
dx=0.5;         %variar para ver estabilidade
dt=0.1;         %variar para ver estabilidade

%vetores
t=0:dt:tf;
x=0:dx:L;

Nt=length(t);
Nx=length(x);
T=zeros(Nx,Nt);

T(:,1)=100; %tempo = inicial (1) -> temperatura = 100
T(1,:)=0;  %colocar o valor inicial da extremidade em 0
T(Nx,:)=0;  %colocar a outra extremidade = 0

D=k/(c*p);  %Constante criada so para simplicar a escrita das equações

%%  Utilizando o método de Crack-Nicolson

