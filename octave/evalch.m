function f=evalch(x0,x1,c,x)

  N = length(c) ;
  t = (x - 0.5*(x1 + x0))/(0.5*(x1 - x0)) ;
  th = acos(t) ;
  T = [] ;

  for i=0:N-1
    T = [T; cos(i*th)] ;
  endfor

  f = c(:)'*T ;
