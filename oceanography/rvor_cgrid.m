% calculates relative vorticity from C-grid model output
% output at interior RHO points
%       [rv,xrv,yrv,zrv] = rvor_cgrid(rgrid,u,v)
% supply rgrid structure with 
%       xu,yu,zu & xv,yv,zv  (all matrices) 
%                   & zw & s_w(vectors)

function [rv,xrv,yrv,zrv] = rvor_cgrid(rgrid,u,v)

    xrv = rgrid.xr(2:end-1,2:end-1,end);
    yrv = rgrid.yr(2:end-1,2:end-1,end);
    zrv = rgrid.zr(2:end-1,2:end-1,:);
    
    gridu.xmat = rgrid.xu; gridu.ymat = rgrid.yu; gridu.zmat = rgrid.zu;
    gridv.xmat = rgrid.xv; gridv.ymat = rgrid.yv; gridv.zmat = rgrid.zv;
    
    gridu.s = rgrid.s_rho; gridv.s = rgrid.s_rho;
    
    vx    = horgrad_cgrid(rgrid,gridv,v,1);
    uy    = horgrad_cgrid(rgrid,gridu,u,2);
    
    rv = avg1(avg1(vx-uy,1),2);