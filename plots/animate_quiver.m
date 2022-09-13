% Makes animation (using pcolor). Assumes last dimension is to be looped by default. 
% Else specify index
%       [] = animate(xax,yax,data,index,pausetime)
%           xax,yax - x,y axes (both optional) - can be empty []
%           data - data to be animated - script squeezes data
%           index - dimension to loop through (optional)
%           pausetime - pause(pausetime) (optional)

function [] = animate_quiver(xax,yax,datax,datay,index,pausetime)

    narg = nargin;
    
    if narg == 2
        datax = squeeze(xax);
        datay = squeeze(yax);
        s    = size(datax);
        xax  = 1:s(1); 
        yax  = 1:s(2);
    end
    
    data = squeeze(datax);
    s = size(datax);
       
    if narg <= 3 || ~exist('index','var')
        stop = s(end);
        index = length(s); % set last dimension to loop by default
    else
        stop = s(index);
    end
    
    if narg ~= 5 || isempty(pausetime)
        pausetime = 0.1;
    end
    
    if narg == 1
        xax = 1:s(1); 
        yax = 1:s(2);
    end
    
    if isempty(xax), xax = 1:s(1); end;
    if isempty(yax), yax = 1:s(2); end;
    
    plotdata = double(squeeze(shiftdim(datax,index)));
    
    datamax = nanmax(plotdata(:));
    datamin = nanmin(plotdata(:));
    
    figure;
    for i=1:stop
        clf
        %ncquiverref(xax,yax,datax(:,:,i)',datay(:,:,i)','',0.5,'col');
        contourf(xax,yax,hypot(datax(:,:,i)',datay(:,:,i)')); shading flat 
        hold on
        quiver(xax,yax,datax(:,:,i)',datay(:,:,i)',1,'Color',[.1 .1 .1]);
        xlim([min(xax) max(xax)]);
        ylim([min(yax) max(yax)]);
        colorbar
        caxis([datamin datamax]);
        title(['index = ' num2str(i)]);
        pause(pausetime); 
    end