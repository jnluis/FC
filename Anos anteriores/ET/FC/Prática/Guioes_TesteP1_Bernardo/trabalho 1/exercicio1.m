%1)
%d)
clc;
clear all;
close all;

h = 0.2;
m=150*10^-3;
g=9.8;
t = 0:0.2:1.5;
N = length(t);
v = zeros(N,1);
z = zeros(N,1);
z(1)= 6;

for k=1:N-1
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k)*h;
end

figure()
subplot(1,2,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,2,2)

plot(t,z);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');
%%
%e)
clc;
clear all;
close all;

h = 0.2;
m=150*10^-3;
g=9.8;
t = 0:0.2:1.5;
N = length(t);
v = zeros(N,1);
z = zeros(N,1);
z(1)= 6;

k=1;
while k < N
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k)*h;
    k = k+1;
end


figure()
subplot(1,2,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,2,2)

plot(t,z);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');
%%
%f) duvidas na expresão analítica
clc;
clear all;
close all;

h = 0.2;
m=150*10^-3;
g=9.8;
t = 0:0.2:1.5;
N = length(t);
v = zeros(N,1);
z = zeros(N,1);
z(1)= 6;

%expressão analítica
vz = v(1) - g*t;
zz = zeros(N,1);
zz(1)= 6;


for k=1:N-1
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k)*h;
    zz(k+1) = zz(k) + vz(k)*h;
end

figure()
subplot(1,2,1)

plot(t,v,t,vz);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,2,2)

plot(t,z,t,zz);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');
%%
%g)o grafico da altura da pedra em função do tempo que já está feito logo
%na primeira alínea xD
clc;
clear all;
close all;

h = 0.2;
m=150*10^-3;
g=9.8;
t = 0:0.2:1.5;
N = length(t);
v = zeros(N,1);
z = zeros(N,1);
z(1)= 6;

for k=1:N-1
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k)*h;
    if(z(k) < 0)
        break;
    end
end

%saber o instante que z=0->interp1(valores de x, valores de y, o valor que 
%queremos)
instante = interp1(t,z,0);
instante = t(instante);

%retirar os valores negativos
t = t(1:7);
v = v(1:7);
z = z(1:7);

%gráficos
figure()
subplot(1,2,1)

plot(t,v);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');

subplot(1,2,2)

plot(t,z);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');
%%
%h) duvidas na expresão analítica
clc;
clear all;
close all;

h = 0.2;
m=150*10^-3;
g=9.8;
t = 0:0.2:1.5;
N = length(t);
v = zeros(N,1);
z = zeros(N,1);
z(1)= 6;

%expressão analítica
vz = v(1) - g*t;
zz(1)= 6;
zz = zz(1) - 0.5*g*t.^2;


for k=1:N-1
    v(k+1) = v(k) - g*h;
    z(k+1) = z(k) + v(k)*h;
end

figure()
subplot(1,2,1)

plot(t,v,t,vz);
title('velocidade');
grid on;
ylabel('v/ms^{-1}');
xlabel('t/s');
legend('numérico','analítico');
subplot(1,2,2)

plot(t,z,t,zz);
title('posição');
grid on;
ylabel('m/m');
xlabel('t/s');
legend('numérico','analítico');

