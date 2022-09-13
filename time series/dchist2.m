function [] = dchist2(x,y,n)

    xi = linspace(min(x(:)),max(x(:)),n);
    yi = linspace(min(y(:)),max(y(:)),n);
    
    xr = interp1(xi,1:numel(xi),x,'nearest');
    yr = interp1(yi,1:numel(yi),y,'nearest');
    
    Z = accumarray([xr yr],1,[n n]);

    surf(Z);
    beautify
    
    figure;
    pcolor(xi,yi,Z);
    beautify;