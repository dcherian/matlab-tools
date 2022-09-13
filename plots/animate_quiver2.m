% Makes animation (using pcolor). Assumes last dimension is to be looped by default. 
% Else specify index
%       [] = animate_quiver2(ax1,data1,ax2,data2,title1,title2,index,pausetime)
%
%           ax* - [xax yax]
%           data% - [datax datay]
%           data - data to be animated - script squeezes data
%           index - dimension to loop through (optional)
%           pausetime - pause(pausetime) (optional)

function [] = animate_quiver2(ax1,data1,ax2,data2,title1, title2, index,pausetime)

    narg = nargin;
    
    if narg <4
        fprintf('\n ERROR: animate_quiver2 requires more arguments.\n\n');
        return
    end
    
    data = squeeze(datax);
    s = size(datax);
       
    if narg <= 3 || ~exist('index','var')
        stop = s(end);
        index = length(s); % set last dimension to loop by default
    else
        stop = s(index);
    end
    
    if narg < 7 || isempty(pausetime)
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
        subplot(1,2,1)
        %ncquiverref(xax,yax,datax(:,:,i)',datay(:,:,i)','',0.5,'col');
        contourf(xax1,yax1,hypot(datax1(:,:,i)',datay1(:,:,i)')); shading flat 
        hold on
        quiver(xax1,yax1,datax1(:,:,i)',datay1(:,:,i)',1,'Color',[.1 .1 .1]);
        xlim([min(xax1) max(xax1)]);
        ylim([min(yax1) max(yax1)]);
        colorbar
        caxis([datamin1 datamax1]);
        title([title1 ' index = ' num2str(i)]);
        
        subplot(1,2,2)
        contourf(xax2,yax2,hypot(datax2(:,:,i)',datay2(:,:,i)')); shading flat 
        hold on
        quiver(xax2,yax2,datax2(:,:,i)',datay2(:,:,i)',1,'Color',[.1 .1 .1]);
        xlim([min(xax2) max(xax2)]);
        ylim([min(yax2) max(yax2)]);
        colorbar
        caxis([datamin2 datamax2]);
        title([title2 ' index = ' num2str(i)]);
        
        pause(pausetime); 
    end