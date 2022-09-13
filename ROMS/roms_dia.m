%function [] = roms_dia(fname)

    fname = 'runs/ocean_dia.nc';
    u_accel  = ncread(fname,'u_accel');
    u_cor    = ncread(fname,'u_cor');
    u_hadv   = ncread(fname,'u_hadv');
    u_prsgrd = ncread(fname,'u_prsgrd');
    
    %%
    ix = 80; iy = 70; iz = 35;
    ut = squeeze(u_accel(ix,iy,iz,:));
    fv = squeeze(u_cor(ix,iy,iz,:));
    uu = squeeze(u_hadv(ix,iy,iz,:));
    px = squeeze(u_prsgrd(ix,iy,iz,:));
    
    bal = ut-fv-uu-px;
    
    t = [1:length(ut)];
    
    plot(t,ut,t,fv,t,uu,t,px,t,bal);
    legend('u_t','fv','nonlinear','p_x','sum');
    