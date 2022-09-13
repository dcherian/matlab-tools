function [] = roms_diagnose_spinup(fname)

    u = double(ncread(fname,'u'));
    v = double(ncread(fname,'v'));
    zeta = double(ncread(fname,'zeta'));
    
    u = avg1(u(:,2:end-1,:,:),1);
    v = avg1(v(2:end-1,:,:,:),2);
    
    [xax,yax,zax,tax,~,~] = dc_roms_var_grid(fname,'rho');
    xax = xax(2:end-1,1,1);
    yax = yax(1,2:end-1,1)';
    
    zax = zax(2:end-1,2:end-1,:);
    
    KE = domain_integrate(0.5 *(u.^2 + v.^2),xax,yax,zax);
    PE = squeeze(trapz(yax,trapz(xax,0.5*9.81*(zeta(2:end-1,2:end-1,:)).^2,1),2));
    
    plot(tax/86400,PE./PE(end) + 1,'b'); hold on
    plot(tax/86400,KE./KE(end),'r');
    xlabel('time (days)');
    legend('PE/PE(end) + 1','KE/KE(end)','Location','SouthEast');