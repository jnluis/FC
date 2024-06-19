%1)
%a)método de euler
clc;
clear all;
close all;

%variáveis
g=9.8;
K=1;
m=1;
w = sqrt(K/m);

%h pode mudar
h=0.01;
tfim=100;

t=0:h:tfim;

N=length(t);

y=zeros(N,1);
y(1)=1;
v=zeros(N,1);
v(1)=0;

%solução numérica
for k=1:N-1
    v(k+1)=v(k) + ( (-K * y(k) ) / m )*h;
    y(k+1)=y(k) + v(k)*h;
end

%gráficos
figure()
subplot(1,3,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,3,2)

plot(t,y);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');

subplot(1,3,3)

plot(y,v);
title('vx em função de x');
grid on;
ylabel('v/ms^{-1}');
xlabel('m/m');

%energia total
% ec = 0.5 * m * v.^2;
% ep = m * y * g;
% em = ec+ep;
et = 0.5 * m * w.^2 * y.^2 + 0.5*m* v.^2;
figure()
plot(t,et);

%%
%b)metodo de euler-cromer
clc;
clear all;
close all;

%variáveis
g=9.8;
K=1;
m=1;
w = sqrt(K/m);
%h pode mudar
h=0.01;
tfim=100;

t=0:h:tfim;

N=length(t);

y=zeros(N,1);
y(1)=1;
v=zeros(N,1);
v(1)=0;

%solução numérica
for k=1:N-1
    v(k+1)=v(k) + ( (-K * y(k) ) / m )*h;
    y(k+1)=y(k) + v(k+1)*h;
end

%gráficos
figure()
subplot(1,3,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,3,2)

plot(t,y);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');

subplot(1,3,3)

plot(y,v);
title('vx em função de x');
grid on;
ylabel('v/ms^{-1}');
xlabel('m/m');

%energia total
% ec = 0.5 * m * v.^2;
% ep = m * y * g;
% em = ec+ep;
et = 0.5 * m * w.^2 * y.^2 + 0.5*m* v.^2;
figure()
plot(t,et);
%%
%c)metodo euler implicito
clc;
clear all;
close all;

%variáveis
g=9.8;
K=1;
m=1;
w = sqrt(K/m);

%h pode mudar
h=0.01;
tfim=100;

t=0:h:tfim;

N=length(t);

y=zeros(N,1);
y(1)=1;
v=zeros(N,1);
v(1)=0;

A =[1 -h; w^2*h 1];
%solução numérica
for k=1:N-1
    %está nos slides teóricos
    b=[y(k) ; v(k)];
    aux = linsolve(A,b);
    y(k+1)= aux(1);
    v(k+1)=aux(2);
end

%gráficos
figure()
subplot(1,3,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,3,2)

plot(t,y);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');

subplot(1,3,3)

plot(y,v);
title('vx em função de x');
grid on;
ylabel('v/ms^{-1}');
xlabel('m/m');

%energia total
% ec = 0.5 * m * v.^2;
% ep = m * y * g;
% em = ec+ep;
et = 0.5 * m * w.^2 * y.^2 + 0.5*m* v.^2;
figure()
plot(t,et);

%%
%d)metodo de crank-nicolson
clc;
clear all;
close all;

%variáveis
g=9.8;
K=1;
m=1;
w = sqrt(K/m);

%h pode mudar
h=0.01;
tfim=100;

t=0:h:tfim;

N=length(t);

y=zeros(N,1);
y(1)=1;
v=zeros(N,1);
v(1)=0;

A =[1 -h*0.5; w^2*h*0.5 1];
%solução numérica
for k=1:N-1
    %está nos slides teóricos
    b=[y(k)+h*0.5*v(k) ; v(k)-w^2*0.5*h*y(k)];
    aux = linsolve(A,b);
    y(k+1)= aux(1);
    v(k+1)=aux(2);
end

%gráficos
figure()
subplot(1,3,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,3,2)

plot(t,y);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');

subplot(1,3,3)

plot(y,v);
title('vx em função de x');
grid on;
ylabel('v/ms^{-1}');
xlabel('m/m');

%energia total
% ec = 0.5 * m * v.^2;
% ep = m * y * g;
% em = ec+ep;
et = 0.5 * m * w.^2 * y.^2 + 0.5*m* v.^2;
figure()
plot(t,et);
%%
%e)
clc;
clear all;
close all;

%variáveis
g=9.8;
K=1;
m=1;
w = sqrt(K/m);

%h pode mudar
h=0.01;
tfim=100;

t=0:h:tfim;

N=length(t);

y=zeros(N,1);
y(1)=1;
v=zeros(N,1);
v(1)=0;

A =[1 -h*0.5; w^2*h*0.5 1];
%solução numérica
for k=1:N-1
    %está nos slides teóricos
    b=[y(k)+h*0.5*v(k) ; v(k)-w^2*0.5*h*y(k)];
    aux = linsolve(A,b);
    y(k+1)= aux(1);
    v(k+1)=aux(2);
end

%solução analítica
y0=y(1);
yx = y0*cos(w.*t);
vx = -y0*w*sin(w.*t);

%para conseguir o periodo e a amplitude
ymax = [];
tmax = [];
imax = 0;

for k = 2:N-1
    if and(y(k+1)-y(k)<=0,y(k)-y(k-1)>=0) % condição de máximo
        imax = imax+1;
        POL = lagr(t(k-1:k+1),y(k-1:k+1));
        tmax(imax) = POL(1);
        ymax(imax) = POL(2);
    end
end

A = mean(ymax);
NI = length(tmax);

PLF = polyfit(1:NI,tmax,1);
figure()
plot(1:NI,tmax)
T = PLF(1);