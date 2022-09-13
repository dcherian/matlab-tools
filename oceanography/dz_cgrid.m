%   Calculates z derivative on a stretched grid
%       [dvdz] = dz_cgrid(grid,var)
%           grid is a structure with s (s for var), s_w, zw

function [dvdz] = dz_cgrid(grid,var)

    dvds = nan(size(var));
    dvds(:,:,2:end-1,:) = avg1(bsxfun(@rdivide,diff(var,1,3),permute(diff(grid.s'),[3 2 1])),3);
    
    % pad on bottom and surface values - IN S
    ds1 = (grid.s(2)-grid.s(1))./(grid.s(3)-grid.s(2));
    ds2 = (grid.s(end)-grid.s(end-1))./(grid.s(end-1) - grid.s(end-2));
    dvds(:,:,1,:)   = dvds(:,:,2,:) - (dvds(:,:,3,:) - dvds(:,:,2,:))*ds1;
    dvds(:,:,end,:) = dvds(:,:,end-1,:) + (dvds(:,:,end-1,:)-dvds(:,:,end-2,:))*ds2;

    % second condition needed if var is w
    if size(grid.zw,1)-1 == size(var,3) || size(grid.zw,1) == size(var,3)
        grid.zw = permute(grid.zw,[3 2 1]);
    end

    dzds = bsxfun(@rdivide,diff(grid.zw,1,3),permute(diff(grid.s_w',1,1),[3 2 1]));

    if size(dzds,1) ~= size(dvds,1)
        dzds = avg1(dzds,1);
    end
    if size(dzds,2) ~= size(dvds,2)
        dzds = avg1(dzds,2);
    end
    if size(dzds,3) ~= size(dvds,3) % if input is w
        dvds = avg1(dvds,3);
    end

    dvdz = bsxfun(@rdivide,dvds,dzds);
WRONNNNNNNNNNNNNNNG!!!!!!!
    % pad on bottom and surface values - IN Z - seems to work better than s
    dz1 = (grid.zmat(:,:,2)-grid.zmat(:,:,1))./(grid.zmat(:,:,3) - grid.zmat(:,:,2));
    dz2 = (grid.zmat(:,:,end)-grid.zmat(:,:,end-1))./(grid.zmat(:,:,end-1) - grid.zmat(:,:,end-2));

    dvdz(:,:,1,:)   = dvdz(:,:,2,:) - ...
                        bsxfun(@times,dvdz(:,:,3,:) - dvdz(:,:,2,:),dz1);
    dvdz(:,:,end,:) = dvdz(:,:,end-1,:) + ...
                        bsxfun(@times,dvdz(:,:,end-1,:)-dvdz(:,:,end-2,:),dz2);

    debug = 0;

    if debug
        range = 1:size(dvdz,3);
        figure
        ax(1) = subplot(141);
        pcolorcen(squeeze(grid.ymat(1,:,range))/1000,squeeze(grid.zmat(1,:,range)),squeeze(1./dzds(1,:,:)));
        shading flat; title('1/dzds = 1/height of grid cell'); colorbar
        ax(2) = subplot(142);
        pcolorcen(squeeze(grid.ymat(1,:,range))/1000,squeeze(grid.zmat(1,:,range)),log(abs(squeeze(dvds(1,:,:)))));
        title('log d/ds'); colorbar;shading flat;
        ax(3) = subplot(143);
        pcolorcen(squeeze(grid.ymat(1,:,range))/1000,squeeze(grid.zmat(1,:,range)),log(abs(squeeze(dvdz(1,:,:)))));
        title('log d/dz = log 1/dzds * d/ds'); colorbar; shading flat
        ax(4) = subplot(144);
        pcolorcen(squeeze(grid.ymat(1,:,range))/1000,squeeze(grid.zmat(1,:,range)),squeeze(zeros(size(grid.xmat(1,:,range)))));
        shading faceted
        linkaxes(ax,'xy');
    end
