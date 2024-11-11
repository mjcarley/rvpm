function [x,w,s]=vtyread(file)

  fid = fopen(file, "r") ;

  np = fscanf(fid, "%d", 1)
  box = fscanf(fid, "%f", 6)
  smax = fscanf(fid, "%f", 1)

  dat = fscanf(fid, "%f") ;

  dat = reshape(dat, 7, np)' ;
  x = dat(:,1:3) ; w = dat(:,4:6) ; s = dat(:,7) ;
  
  fclose(fid) ;
