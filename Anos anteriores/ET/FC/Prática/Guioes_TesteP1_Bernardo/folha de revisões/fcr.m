function F = fcr(xv,xold,vold,const)
    F(1) = xv(1)-xold-(vold+xv(2))*const(1);
    F(2) = xv(2)-vold+const(2)*(xold+xv(1)+const(3)*(xold^3+xv(1)^3));
end