function derivadas = f(t,solucao,K,m)
    derivadas = zeros(2,1);
    x = solucao(1);
    v = solucao(2);
    derivadas(1) = v;
    derivadas(2) = -K/m*x;
end