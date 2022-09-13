% Returns Beckmann & haidvogel; Haney numbers in that order
%           [r_x0,r_x1,rx0mat,rx1mat] = stiffness(h,zrmat)

function [r_x0,r_x1,rx0mat,rx1mat] = stiffness(h,zrmat)

    xflag = 0; yflag = 0;
    if bsxfun(@minus,h,h(:,1)) == zeros(size(h)), xflag = 1; end
    if bsxfun(@minus,h,h(1,:)) == zeros(size(h)), yflag = 1; end
    
    % Grid metrics - From Utility/stiffness.F
    % beckmann & haidvogel (1993) number 
    if xflag
        r_x0_X = abs(diff(h,1,1)./avg1(h,1))/2; 
    else
        r_x0_X = 0;
    end
    if yflag
        r_x0_Y = abs(diff(h,1,2)./avg1(h,2))/2;
    else
        r_x0_Y = 0;
    end
    r_x0 = max([max(r_x0_X(:)) max(r_x0_Y(:))]);

    rx0mat = max(avg1(r_x0_X,2),avg1(r_x0_Y,1));
    
    % haney (1991) number
    r_x1_x = nan(size(zrmat));
    r_x1_y = nan(size(zrmat));

    for k=2:size(zrmat,3)
        if xflag
          for i=2:size(zrmat,1)
            r_x1_x(i,:,k) = (zrmat(i,:,k) - zrmat(i-1,:,k) + zrmat(i,:,k-1) - zrmat(i-1,:,k-1)) ...
                            ./ (zrmat(i,:,k) + zrmat(i-1,:,k) - zrmat(i,:,k-1) - zrmat(i-1,:,k-1));
          end
        else
            r_x1_x = 0;
        end
        if yflag
          for j=2:size(zrmat,2)

              r_x1_y(:,j,k) = (zrmat(:,j,k) - zrmat(:,j-1,k) + zrmat(:,j,k-1) - zrmat(:,j-1,k-1)) ...
                            ./ (zrmat(:,j,k) + zrmat(:,j-1,k) - zrmat(:,j,k-1) - zrmat(:,j-1,k-1));
          end
        else
            r_x1_y = 0;
        end
    end
    r_x1 = nanmax(abs( [r_x1_x(:);r_x1_y(:)] ));
    rx1mat = max(abs(r_x1_x),abs(r_x1_y));