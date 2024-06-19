% Parte 2
%*************************************************************************
%
% NOME 1: Martim Gil
% MEC  1: 102901
% Turma : PL6
% 
%*************************************************************************
%
% NOME 2: João Luís
% MEC  2: 107403
% Turma : PL6
% 
%*************************************************************************
%
% NOME 3: João Marques
% MEC  3: 108072
% Turma : PL6
%
%*************************************************************************
%
% NOME 4: Délcio Amorim
% MEC  4: 109680
% Turma : PL6
%
%*****************************Alínea e)***********************************
clc;
clear all;
close all;

colors = ['b','g','r'];
nn=[1 3 0.5];
for is=1:3

N=1024;
L=80;
dx=L/N;
alfa=6;
n=nn(is);
C=10;
x=-(N-1)/2*dx:dx:(N-1)/2*dx;

ti=0;
tf=1;
dt=1.5e-5;
t=ti:dt:tf;
Nt=length(t);

% Vetor das frequencias angulares
dw=2*pi/(N*dx);
wmax=(N/2-1);
wmin=(-N/2);
w=[0:wmax wmin:-1]*dw;

qx=zeros(Nt,N);
qx(1,:)=((C/2)*((sech((sqrt(C)*n.*x)/2)).^(2))).^(1/n);
%qx(1,:)=((C/2)*(sech(sqrt(C))*n.*x/2).^(2)).^(1/n);
t1=(1i.*w);
t3=(1i.*w).^3;

for r = 1:Nt-1
    
    q=qx(r,:);
    
    r1 = (-ifft(t3.*fft(q)) - (n+1)*(n+2).*q.^(n).*ifft(t1.*fft(q)));
    v = q + r1*dt/2;
    r2 = (-ifft(t3.*fft(v)) - (n+1)*(n+2).*v.^(n).*ifft(t1.*fft(v)));
    v2 = q + r2*dt/2;
    r3 = (-ifft(t3.*fft(v2)) - (n+1)*(n+2).*v2.^(n).*ifft(t1.*fft(v2)));
    v3 = q + r3*dt;
    r4 = (-ifft(t3.*fft(v3)) - (n+1)*(n+2).*v3.^(n).*ifft(t1.*fft(v3)));
   
    qx(r+1,:)= qx(r,:) + 1/6*(r1 + 2*r2 + 2*r3 + r4)*dt;
             
end
a=(qx(end,:).^(2))';
% Encontrar e imprimir os picos
velocidade = [];
posicao = [];
amplitude = [];
for i = 2:N-1
    if a(i) > a(i-1) && a(i) > a(i+1)
        if (a(i)>=1)
        posicao = [posicao, x(i)];
        amplitude = [amplitude, a(i)];
        velocidade = [velocidade,(x(i)-0)/tf];
        end
    end
end
fprintf('para n = %.1f\n', n);
for i = 1:length(posicao)
    fprintf('Pico %d: posição = %.2f, Amplitude = %.2f , velocidade = %.2f\n', i, posicao(i), amplitude(i),velocidade(i));
end

plot(x,q(end,:),colors(is), 'DisplayName', ['n = ' num2str(n)]);
hold on
end
title('Perfis finais em função de x')
xlabel('x')
ylabel('q')
legend;
hold off
