%%  Ex1.c

clear,clc,close all

%% CONSTANTES

h=0.01;
tf=10;
t=0:h:tf;

k=16;
m=1;

x=[];
v=[];

x(1)=1;
y(1)=10;

%% MÉTODO DE RUNGE-KUTTA

