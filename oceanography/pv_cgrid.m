% Calculate PV on C-grid. Result at interior RHO points
%        [pv,xpv,ypv,zpv] = pv_cgrid(rgrid,u,v,rho,f,rho0)
% supply rgrid structure with 
%       xu,yu,zu & xv,yv,zv & xr,yr,zr (all matrices) 
%                   & zw & s_w(vectors), s_rho

function [pv,xpv,ypv,zpv,rvor] = pv_cgrid(rgrid,u,v,rho,f,rho0)

    xpv = avg1(rgrid.xr(2:end-1,2:end-1,:), 3);
    ypv = avg1(rgrid.yr(2:end-1,2:end-1,:), 3);
    zpv = avg1(rgrid.zr(2:end-1,2:end-1,:), 3);
    
    gridu.xmat = rgrid.xu; gridu.ymat = rgrid.yu; gridu.zmat = rgrid.zu;
    gridv.xmat = rgrid.xv; gridv.ymat = rgrid.yv; gridv.zmat = rgrid.zv;
    gridr.xmat = rgrid.xr; gridr.ymat = rgrid.yr; gridr.zmat = rgrid.zr;
    
    gridu.s = rgrid.s_rho; gridv.s = rgrid.s_rho; gridr.s = rgrid.s_rho;
    gridu.zw = rgrid.zw; gridv.zw = rgrid.zw; gridr.zw = rgrid.zw; 
    gridu.s_w = rgrid.s_w; gridv.s_w = rgrid.s_w; gridr.s_w = rgrid.s_w;
    
    s = size(rho);
    if length(s) == 3

        s(4) = 1;
    end
    
    if length(f) == 1 % constant f0
        f   = repmat(f,[1 1 s(3) s(4)]);
    end
    
    % average to land on vx - uy points
    f = avg1(avg1(avg1(f,2),2),1);
        
    if ~exist('rho0','var');
        fprintf('\n\n Using rho0 = 1027. \n\n');
        rho0 = 1027;
    end
    
    % calculate gradients
    vx    = diff_cgrid(gridv,v,1); 
    %vy    = diff_cgrid(gridv,v,2); 
    vz    = diff_cgrid(gridv,v,3);

    %ux    = diff_cgrid(gridu,u,1); 
    uy    = diff_cgrid(gridu,u,2); 
    uz    = diff_cgrid(gridu,u,3);

    tx    = diff_cgrid(gridr,rho,1); 
    ty    = diff_cgrid(gridr,rho,2); 
    tz    = diff_cgrid(gridr,rho,3);
    
    if size(tz,3) > size(vx,3)
        tz = avg1(tz,3);
        uz = avg1(uz,3);
        vz = avg1(vz,3);
    end
    
    rvor  = vx - uy;

    % PV calculated at interior rho points
                                % f + vx - uy                      (rho)_z
    pv = -1* double((avg1(bsxfun(@plus,avg1(rvor,2),f), 1).*tz(2:end-1,2:end-1,:,:) ...
                   - avg1(vz(2:end-1,:,:,:),2).*avg1(tx(:,2:end-1,:,:),1) ... % vz * (rho)_x
                   + avg1(uz(:,2:end-1,:,:),1).*avg1(ty(2:end-1,:,:,:),2))./rho0);%avgz(lambda(2:end-1,2:end-1,:,:))); % uz*(rho)_y
               
    debug = 0;
    if debug
        pv1 = -avgx(avgz(bsxfun(@plus,avgy(vx - uy),f)))  .*  tz(2:end-1,2:end-1,:,:);
        pv2 = avgy(vz(2:end-1,:,:,:)).*avgz(avgx(tx(:,2:end-1,:,:)));
        pv3 = avgx(uz(:,2:end-1,:,:)).*avgz(avgy(ty(2:end-1,:,:,:)));
        
        tind = 1;
        yind = 3;
        
        figure;
        contourf(squeeze(xpv(:,yind,:)),zpv,squeeze(pv1(:,yind,:,tind))',20);colorbar
        title('(f + v_x -u_y)\rho_z');
        figure;
        contourf(xpv,zpv,squeeze(pv2(:,yind,:,tind))',20);colorbar
        title('v_z \rho_x');
        figure;
        contourf(xpv,zpv,squeeze(pv3(:,yind,:,tind))',20);colorbar
        title('u_x \rho_y');
        figure;
        contourf(xpv,zpv,squeeze(pv(:,yind,:,tind))',20);colorbar
        title('Full PV');
        colormap(hsv);
        pause;
    end
        
% following for *flat bottom* only
%     vx    = bsxfun(@rdivide,diff(v,1,1),diff(grid.xv)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
%     vy    = bsxfun(@rdivide,diff(v,1,2),diff(grid.yv')); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
%     ux    = bsxfun(@rdivide,diff(u,1,1),diff(grid.xu)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
%     uy    = bsxfun(@rdivide,diff(u,1,2),diff(grid.yu')); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
%     tx    = bsxfun(@rdivide,diff(rho,1,1),diff(grid.xr)); %diff(v,1,1)./repmat(diff(grid.x_v',1,1),[1 1 s(3) s(4)]);
%     ty    = bsxfun(@rdivide,diff(rho,1,2),diff(grid.yr')); %diff(v,1,2)./repmat(diff(grid.y_v',1,2),[1 1 s(3) s(4)]);
    
