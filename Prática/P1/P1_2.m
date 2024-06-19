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
for k=1:N-1
    dV(k+1) = dV(k)+ (-a*V(k)) *h;
    V(k+1) = V(k) + dV(k) * h; % EULER
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
legend("numérico","analítico");
xlabel("Tempo (s)")
ylabel("Tensão (V)")

figure(2);
plot(t,Q,t, Qt);
legend("numérico","analítico");
xlabel("Tempo (s)")
ylabel("Carga (Coloumb)")

figure(3);
plot(t,IL,t, Ilt); % a corrente está desfasada da tensão e da carga por pi/2
legend("numérico","analítico");
xlabel("Tempo (s)")
ylabel("Corrente (Ampére)")

%% Alinea g

% Valor teórico do período é dado por T= (2*pi)/ W
Tteorico = (2*pi)/ W

Maxlocals = islocalmax(V);
tMaxlocals = t(Maxlocals); % Instantes dos máximos locais

II = find(islocalmax(V)>0); %% forma alternativa ao Maxlocals, dá igual
TT=t(II);

nI= length(II);

%% PERIODO com interpolaçao
for ii=1:nI
    j=II(ii);
    V00=V(j);  
    Tm(ii)=interp1(V(j-1:j+1),t(j-1:j+1),V00,'linear');
end

% Este trecho utiliza um loop for para iterar sobre os índices de máximos locais (II).
% Para cada índice j, ele utiliza a função interp1 para realizar interpolação linear nos valores de t usando os valores vizinhos de V.
% O resultado é armazenado no vetor Tm, que representa os tempos interpolados associados aos máximos locais.

for ii=2:nI
    jj=ii-1;  
    Texp2(jj)=Tm(ii)-Tm(ii-1);
end

% Este trecho acima utiliza outro loop for para calcular as diferenças
% de tempo entre os tempos interpolados consecutivos e armazenar essas diferenças no vetor Texp2.

TE2=mean(Texp2) %% PERIODO da solução numérica