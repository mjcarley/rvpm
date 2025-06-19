function [K,dK,k,dk,g,dg,Ks]=gskernel(x, y, s)

  r = [x(:,1) - y(:,1) x(:,2) - y(:,2) x(:,3) - y(:,3)] ;

  R2 = sum(r.^2,2) ;
  R = sqrt(R2) ;
  R3 = R.*R2 ;
  
  g = erf(R/s) - 2/sqrt(pi)*R/s.*exp(-R2/s^2) ;

  k = -[r(:,1) r(:,2) r(:,3)]./[R3 R3 R3]/4/pi ;
  
  K = -[r(:,1).*g r(:,2).*g r(:,3).*g]./[R3 R3 R3]/4/pi ;

  R5 = R3.*R2 ;

  dk = [] ;
  for i=1:3
    for j=1:3
      dk = [dk -3*r(:,i).*r(:,j)./R5] ;
    endfor
  endfor

  dk(:,1) += 1./R3 ;
  dk(:,5) += 1./R3 ;
  dk(:,9) += 1./R3 ;

  E = exp(-R2/s^2) ;
  dg = [r(:,1).*R/s^3*4/sqrt(pi).*E r(:,2).*R/s^3*4/sqrt(pi).*E ...
	r(:,3).*R/s^3*4/sqrt(pi).*E] ;

  dK = [] ;
  for i=1:3
    for j=1:3
    dK = [dK -dk(:,3*i+j-3).*g/4/pi+k(:,j).*dg(:,i)] ;
    endfor
  endfor

  i = find(R == 0) ;
  K(i,:) = 0 ;

  dK(i,:) = 0 ;
  dK(i, [1 5 9]) = -1/pi/sqrt(pi)/3/s^3 ;

  ## series expansion
  Ks = -[r/s^3].*[E E E]/3/pi^1.5 - [4*r/s^5/15].*[R.^2.*E R.^2.*E R.^2.*E]/pi^1.5 ;
