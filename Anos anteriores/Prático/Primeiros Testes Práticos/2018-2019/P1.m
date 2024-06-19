%% Problema 1.2 - Oscilador harmónico simples

clc;
close all;
clear all;

L = 0.1; % indutância
C = 1e-3; % capacidade
epsilon = 5; % V
R= 10; % Ohms

a = 1 / (R*C) + R/L;
b = 2 / (L *C);
c=epsilon / (L*C);

t0 = 0;
tf= 0.5;
h= 0.001;
t = t0:h:tf;

N= numel(t); % igual ao length, porque não há matrizes mas apenas array's
V = zeros(1,N); % VC
dV = zeros(1,N);

dV(1) =0; % Nos dados

fdV = @(V, dV) -a*dV - b*V + c; % função anónima para dV/dt
fV = @(dV) dV; % função anónima para V

for k=2:N-1
    % EULER
    %dV(k+1) = dV(k)+ (-a*dV(k) - b*(V(k) )+c) *h;
    %V(k+1) = V(k) + dV(k) * h; 

    %Implícito com LinSolve
    % Definindo a matriz A e o vetor B para o sistema A*X = B
    % onde X é o vetor [dV(k+1); V(k+1)]
    % e B depende de dV(k), V(k) e possivelmente de valores constantes
    
    % A = [1, -h; b*h, 1+a*h];
    % B = [V(k); dV(k) + h*c];
    % 
    % % Usando linsolve para resolver A*X = B
    % X = linsolve(A, B);
    % 
    % dV(k+1) = X(1);
    % V(k+1) = X(2);

    % Runge-Kutta
    r1V = h*fV(dV(k));
    r1dV = h*fdV(V(k), dV(k));
    
    r2V = h*fV(dV(k) + r1dV/2);
    r2dV = h*fdV(V(k) + r1V/2, dV(k) + r1dV/2);
    
    V(k+1) = V(k) + r2V;
    dV(k+1) = dV(k) + r2dV;

    if and(abs(V(k) - V(k-1)) < 1e-6 , abs(V(k) - (epsilon/2)) < 1e-6)
        break
    end
end

t = t(1:k+1); % descarta zeros
V = V(1:k+1);
dV = dV(1:k+1);

Q = V *C; % Carga
IL = dV * C; % Corrente

W= sqrt(1/(L*C));
Q0 = 5*C;
ILMax= W*Q0; 

% Solução analitica
% Vt = V0* cos(W*t);
% Qt= Q0* cos(W*t);
% Ilt= -W*Q0* sin(W*t);

figure(1);
plot(t,V, 'r');
xlabel("Tempo (s)")
ylabel("Vc")
N=k;
% Localiza máximos
imax = 0;
for k= 2:N-2
    if and(V(k+1)-V(k) <=0 , V(k) - V(k-1) >= 0 )
        imax = imax+1;
        aux = lagr(t(k-1:k+1), V(k-1:k+1));
        tmax(imax) = aux(1);
        Vmax(imax) = aux(2);
    end
end

% Encontra o maior máximo local
if imax > 0
    [VmaxGlobal, idxMaxGlobal] = max(Vmax);
    tmaxGlobal = tmax(idxMaxGlobal);
    fprintf('Maior valor máximo de Vc método numérico: %f V\n', VmaxGlobal);
    fprintf('Instante de ocorrência do maior valor máximo: %f s\n', tmaxGlobal);
else
    fprintf('Nenhum máximo encontrado.\n');
end