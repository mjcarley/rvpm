z = [] ;
t = [] ;
dt = 0.0125 ;

##dat=load(["points.dat"]) ;
##x=dat(:,1:3) ; w=dat(:,4:6) ; s = dat(:,7) ;
##z = [z; mean(x(:,3))] ;

for i = 0:10:190
  dat=load(["points-" int2str(i) ".dat"]) ;
  x=dat(:,1:3) ; w=dat(:,4:6) ; s = dat(:,7) ;

  z = [z; mean(x(:,3))] ;
  t = [t; (i+1)*dt] ;
  
endfor
