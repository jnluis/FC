clc;
clear all;
close all;

T= 1000;
omega = 7741.8;% frequencua
L= 0.35; %m
vsom=345; % velocidade do som

xf=L;
h=0.0001;
x=0:h:xf;
N=numel(x);

y=zeros(N,1);
Dy=zeros(N,1);

y(1)=5*10^-6;
Dy(1)=0;

for k=1:N-1
    Dy(k+1)=Dy(k)+ (((-omega^2) / vsom) * y(k)) *h;
    y(k+1) = y(k)+Dy(k+1)*h;
end

plot(x,Dy);
xlabel('x'); ylabel('y');

n1=1;
n3=3;
n5=5;
freqpropr1 = ((n1* pi) / (2*L)) * vsom
freqpropr3 = ((n3* pi) / (2*L)) * vsom
freqpropr5 = ((n5* pi) / (2*L)) * vsom