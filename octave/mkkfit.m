## generate a Chebyshev interpolant for Gaussian kernel
tol = 1e-13 ;

## cut off for kernel evaluation (beyond this point equal to 1 to
## machine precision)
xmax = 6.3 ;

## order of Chebyshev polynomials to use, roots and values
N = 10 ;
t = cos(((N-1:-1:0)+1/2)*pi/N) ;
th = acos(t) ;
T = [] ;
for i=0:N-1
  T = [T; 2*cos(i*th)] ;
endfor
T(1,:) *= 0.5 ;

## interpolation intervals
xi = linspace(0,xmax,2) ;

## test points
tc = linspace(-1,1,65) ;

C = [] ;
i = 0 ;
while ( i < length(xi)-1)
  ##for i=1:length(xi)-1
  i ++ ;
  x = 0.5*(xi(i+1) + xi(i)) + 0.5*(xi(i+1) - xi(i))*t ;
  g = gfunc(x) ;
  ##g = x.^3 ;
  c = T*g'/N ;

  ## check for convergence
  xc = 0.5*(xi(i+1) + xi(i)) + 0.5*(xi(i+1) - xi(i))*tc ;
  g = gfunc(xc) ;
  f = evalch(xi(i), xi(i+1), c, xc) ;
  if ( max(abs(f-g)) > tol )
    xi = [xi(1:i) 0.5*(xi(i+1)+xi(i)) xi(i+1:end)] ;
    i -- ;
  else
    C = [C c] ;
  endif
endwhile

x = linspace(0, xmax-1e-3, 10000) ;
f = 0*x ;
for i=1:length(x)
  ii = find(xi(1:end-1)<=x(i) & xi(2:end)>x(i)) ;
  f(i) = evalch(xi(ii), xi(ii+1), C(:,ii), x(i)) ;
endfor

g = gfunc(x) ;
if ( max(abs(f-g)) > tol )
  error("tolerance not achieved") ;
endif

fid = fopen("gkernel-data.c", "w") ;

fprintf(fid, "#include <stdlib.h>\n\n") ;
fprintf(fid, "#include <glib.h>\n\n") ;

fprintf(fid, "#define RVPM_GAUSS_FIT_ORDER %d\n\n", N) ;
fprintf(fid, "#define RVPM_GAUSS_FIT_CUTOFF %1.16e\n\n", xmax) ;
fprintf(fid, "#define RVPM_GAUSS_FIT_INTERVAL_NUMBER %d\n\n", length(xi)-1) ;

fprintf(fid, "gdouble RVPM_GAUSS_FIT_COEFFICIENTS[] = {\n") ;
dat = C ;
fprintf(fid, "  %1.16e,\n", dat(1:end-1)) ;
fprintf(fid, "  %1.16e} ;\n\n", dat(end)) ;

fprintf(fid, "gdouble RVPM_GAUSS_FIT_INTERVALS[] = {\n") ;
dat = xi' ;
fprintf(fid, "  %1.16e,\n", dat(1:end-1)) ;
fprintf(fid, "  %1.16e} ;\n", dat(end)) ;

fclose(fid) ;
