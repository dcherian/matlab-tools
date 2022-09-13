% Draws horizontal or vertical line and adds tick mark
% called by linex and liney
function [handles, txthandles] = dcline(ax,x,label,color)

    if isempty(x)
        handles = [];
        txthandles = [];
        return;
    end

    if isempty(label), txthandles = []; end

    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    if size(x,2) == 1, x = x'; end
    
    figure(hFig);
    hold on;
    
    if ~exist('label','var'), label = [];num2str(x'); end
    if ~exist('color','var') || isempty(color), color = [1 1 1]*0.5; end
    
    if length(x) ~= size(label,1)
        label = repmat(label,[length(x) 1]);
    end
    
    if ~iscell(label) && ~isempty(label)
        label = cellstr(label);
    end
    
    if ax == 'x'
        yax  = get(hAxis,'YLim');
        tickstr = 'XTick';
    else
        xax  = get(hAxis,'XLim');
        tickstr = 'YTick';
    end
    
    tick = get(hAxis,tickstr);

    txtfactor = 0.15;

    fs = get(groot, 'DefaultTextFontSize') - 2;

    for i=1:length(x)
        if ax == 'x'
            handles(i) = plot([x(i) x(i)],yax,'-','LineWidth',1,'Color',color);
            if ~isempty(label)
                txthandles(i) = text(double(x(i)),double(yax(1) + txtfactor * diff(yax)), ...
                                     ['  ' label{i}], ...
                                     'Rotation',90,'VerticalAlignment', ...
                                     'Bottom','FontSize',fs,'Color',color);
            end
        else
            handles(i) = plot(xax,[x(i) x(i)],'-','LineWidth',1,'Color',color);
            if ~isempty(label)
                txthandles(i) = text(double(xax(1) + txtfactor*diff(xax)), double(x(i)), ...
                                     ['  ' label{i}], ...
                                     'Rotation',0,'VerticalAlignment','Bottom', ...
                                     'FontSize',fs,'Color',color);
            end
        end

        handles(i).Tag = 'dcline'; % tag so that I can find them later
        handles(i).HandleVisibility = 'off'; % prevent legend from seeing it
    end
    drawnow;
    if length(handles) == 1
        handles = handles(1);
        if ~isempty(label)
            txthandles = txthandles(1);
        end
    end
    % add extra axis ticks
    set(hAxis,tickstr,sort(unique(str2num(sprintf('%.2e ', [tick x])))));
