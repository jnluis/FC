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
%*****************************Alínea d)***********************************
clc;
clear all;
close all;

N=1024;
L=160;
dx=L/N;
alfa=6;

x=-(N-1)/2*dx:dx:(N-1)/2*dx;
%x1=x+60;
x1=x;
ti=-4;
tf=8;
dt=1.5e-5;
t=ti:dt:tf;
Nt=length(t);

% Vetor das frequencias angulares
dw=2*pi/(N*dx);
wmax=(N/2-1);
wmin=(-N/2);
w=[0:wmax wmin:-1]*dw;

qx=zeros(Nt,N);
%qx(1,:)=-12.*sech(x1).^(2);
num=3 + 4.*cosh(2.*x1 + 24*ti) + cosh(4.*x1);
dem=(3.*cosh(x1 - 12*ti) + cosh(3.*x1 + 12*ti)).^(2);
cond_inicial=-12.*num./dem;
qx(1,:)=cond_inicial;
t1=(1i.*w);
t3=(1i.*w).^3;

%%
for n = 1:Nt-1
    
    q=qx(n,:);
    
    r1 = (-ifft(t3.*fft(q)) + alfa.*q.*ifft(t1.*fft(q)));
    v = q + r1*dt/2;
    r2 = (-ifft(t3.*fft(v)) + alfa.*v.*ifft(t1.*fft(v)));
    v2 = q + r2*dt/2;
    r3 = (-ifft(t3.*fft(v2)) + alfa.*v2.*ifft(t1.*fft(v2)));
    v3 = q + r3*dt;
    r4 = (-ifft(t3.*fft(v3)) + alfa.*v3.*ifft(t1.*fft(v3)));
   

    qx(n+1,:)= qx(n,:) + 1/6*(r1 + 2*r2 + 2*r3 + r4)*dt;
       
       
end

%%

a=(qx(end,:).^(2))';
figure();
plot(x,qx(end,:).^(2))
title('Perfil final de q')
xlabel('x')
ylabel('q')
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
fprintf('para alpha = %.1f\n', alfa);
for i = 1:length(posicao)
    fprintf('Pico %d: posição = %.2f, Amplitude = %.2f , velocidade = %.2f\n', i, posicao(i), amplitude(i),velocidade(i));
end
figure()
plot(velocidade,amplitude,'DisplayName', ['\alpha = ' num2str(alfa)]);
title('Amplitude em função da velocidade')
xlabel('velocidade')
ylabel('Amplitude')
legend;

%%
figure()
intensidade=abs(qx).^2;
mesh(x1,t(1:500:end),intensidade(1:500:end,:))
axis([-Inf, Inf, -Inf,Inf,-Inf,Inf])
set(gca,'YDir','reverse');
title('Evolução de q(t, x) usando a equação KdV')
xlabel('x');
ylabel('t')
zlabel('Intensidade')

