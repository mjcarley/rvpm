function K=kernel(x, y)

  r = x - y ;

  R = sqrt(r*r') ;
  if ( R < 1e-6 ) K = [0 0 0] ; return ; endif
  K = -r/R^3/4/pi ;
  
