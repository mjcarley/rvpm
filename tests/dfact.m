function d=dfact(n)

if ( n == -1 ) d = 1 ; return ; end
if ( n == 0 ) d = 1 ; return ; end

if ( 2*floor(n/2) == n ) 
   d = prod(2:2:n) ;
else
  d = prod(1:2:n) ;
end
