function derivadas = f3_3(t, solucao, m,K)
w= sqrt(K/m);

derivadas = zeros(2,1);
x=solucao(1);
v=solucao(2);

derivadas(1)= v;
derivadas(2)= -w^2*x;