y = [0.3 0.1 0.7] ;

x = [linspace(-0.125,0.125,65)' 0*linspace(-0.5,0.5,65)' 0*linspace(-0.5,0.5,65)'] ;

x = [x(:,1)+y(1) x(:,2)+y(2) x(:,3)+y(3)] ;

s = 0.1 ;

[K,dK,k,dk,g,dg,Ks] = gskernel(x, y, s) ;

