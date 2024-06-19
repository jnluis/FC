%% Problema 1.2 - Oscilador harmónico simples

clc;
close all;
clear all;

L = 0.25; % indutância
C = 1e-3; % capacidade
V0 = 5;

t0 = 0;
tf= 0.5;
h= 0.001;
t = t0:h:tf;

a = 1/(L*C);
N= numel(t); % igual ao length, porque não há matrizes mas apenas array's
V = zeros(1,N); % VC
dV = zeros(1,N);

V(1) =V0;
dV(1) =0; % como é uma bobina, a corrente não pode ser descontinuada.
        % Assim, I(t)=0 => derivada de V_C é 0
Uc = zeros(1,N);
Um = zeros(1,N);

for k=1:N-1
    dV(k+1) = dV(k)+ (-a*V(k)) *h;
    V(k+1) = V(k) + dV(k+1) * h; % EULER-CROMER
end

t = t(1:k+1); % descarta zeros
V = V(1:k+1);
dV = dV(1:k+1);

Q = V *C; % Carga
IL = dV * C; % Corrente

W= sqrt(1/(L*C))
Q0 = 5*C;
ILMax= W*Q0 

% Solução analitica
Vt = V0* cos(W*t);
Qt= Q0* cos(W*t);
Ilt= -W*Q0* sin(W*t);

figure(1);
plot(t,V,t,Vt, 'r');
legend("Euler-Cromer","analítico");
xlabel("Tempo (s)")
ylabel("Tensão (V)")

% Energia Elétrica e magnética numérica
Uc = 0.5 * (Q.^2) ./ C;
Um = 0.5 * L * (IL.^2);

Em = Uc + Um;

figure(2);
plot(t,Uc,"r", t, Um, "k");
xlabel("Tempo (s)")
ylabel("Energia Elétrica (J) e Energia Magnética (J)")

% Energia Total Analítica
EmA = 0.5 * (Q0.^2) / C;
EmA = ones(size(Em)) * EmA; % Para aparecer como reta no gráfico em vez de um ponto
figure(3);
plot(t,Em,"r",t, EmA,"b"); 
legend("Energia total (valor numérico)","Energia total (Valor esperado)");
xlabel("Tempo (s)")
ylabel("Energia Total (J)")