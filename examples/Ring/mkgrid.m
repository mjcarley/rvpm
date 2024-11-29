dx = dz = 0.01 ;

xmin = 0.6 ; xmax = 1.2 ; nx = floor((xmax - xmin)/dx) ;
y = 0 ;
zmin = -0.25 ; zmax = 0.5 ; nz = floor((zmax - zmin)/dz) ;

[x,z] = meshgrid(linspace(xmin, xmax, nx), linspace(zmin, zmax, nz)) ;

dat = [x(:) y+0*x(:) z(:)] ;

fid = fopen("grid.dat", "w") ;

fprintf(fid, "%f %f %f\n", dat') ;

fclose(fid) ;

## rvpm-solve -v grid.dat < ring.dat > w.dat
## wy=reshape(dat(:,5), nz, nx)  ;
## contour(z, x, wy, "k")
