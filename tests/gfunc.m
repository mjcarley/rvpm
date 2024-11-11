r = linspace(0, 4, 1000) ;

g = erf(r) - 2/sqrt(pi)*r.*exp(-r.^2) ;
g1 = 1-erfc(r) - 2/sqrt(pi)*r.*exp(-r.^2) ;

G = 0 ;

for k=1:16
  G += 2^k*r.^(2*k)/dfact(2*k+1) ;
endfor

G .*= 2/sqrt(pi)*r.*exp(-r.^2) ;



