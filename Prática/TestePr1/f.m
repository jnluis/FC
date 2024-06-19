function derivadas = f(t, solucao, m,K)
g= 9.81;
d = 0.025;
D= 0.35;

derivadas = zeros(2,1);
x=solucao(1);
v=solucao(2);

derivadas(1)= v;
derivadas(2)= -sqrt(2*g)*((d/D).^2) * sqrt(x);

