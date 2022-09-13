% Makes animation (default using contourf). Assumes last dimension is to be looped by default. 
% Else specify index. Allows browsing.
%       [] = animate(xax,yax,data,labels,commands,index,pausetime)
%           xax,yax - x,y axes (both optional) - can be empty []
%           data - data to be animated - script squeezes data
%
%           labels - Structure with following fields for labeling plot (optional - can be [])
%               > title
%               > xax
%               > yax
%               > revz - reverse yDir?
%               > time - where to pull time info for titles
%               > tunits - units for time vector labels.time
%               > tmax  - maximum value of time vector : needed when only
%                         part of data is being plotted at a time. set by
%                         default to max(labels.time) if not defined
%               > stride - which stride are we on when called by mod_movie? (default: 0)
%               > dt - time step interval from mod_movie (default: 1)
%               > t0 - inital timestep from mod_movie
%
%           commands - custom commands to execute after plot or one of below. (in string, optional)
%                    - separate commands using ; 
%                    - Built-in options are: 
%                          > nocaxis - non-constant colorbar
%                          > pause   - start paused
%                          > fancy_cmap - LAB space colormap
%                          > pcolor  - use pcolor instead of contourf
%                          > imagesc - use imagesc(nan) instead of contourf. imagescnan is tried first
%                          > contour - use contour instead of contourf
%
%           index - dimension to loop through (optional)
%           pausetime - pause(pausetime) (optional)
%
% USAGE:
%       animate(data)
%       animate(data,commands)
%       animate(xax,yax,data)
%       animate([],[],data,...
%       animate(xax,yax,data,commands)
%               eg: animate(xax,yax,data,'pcolor;nocaxis')
%
% BROWSE:
%       - first space pauses, second space resumes, remaining spaces *play*
%       - arrowkeys *pause* and navigate always
%       - Esc to quit

function [] = animate(xax,yax,data,labels,commands,index,pausetime)

    %% figure out input
    narg = nargin;
    
    if strcmp(get(gcf,'currentkey'),'escape')
        warning('Previous ESC detected. Opening new figure.');
        figure;
    end
    
    switch narg
        case 1,
            data = squeeze(xax);
            s    = size(data);
            xax  = 1:s(1); 
            yax  = 1:s(2);
            
        case 2,
            data = xax;
            commands = yax;
            xax = [];
            yax = [];            
            
        case 4,
            if ischar(class(labels))
               commands = labels;
               labels = [];
            end
    end
    
    data = squeeze(data);
    s = size(data);
    if length(s) == 2, s(3) = 1; end
       
    if narg <= 5 || ~exist('index','var')
        stop = size(data,3);
        index = length(s); % set last dimension to loop by default
    else
        stop = s(index);
    end
    
    if narg ~= 7 || ~exist('pausetime','var')
        pausetime = 0.2;
    end
    
    if ~exist('labels','var') || isempty(labels)
        labels.title = '';
        labels.xax   = '';
        labels.yax   = '';
        labels.revz  = 0;
        labels.time  = [];
        labels.tmax  = size(data,3);
        labels.stride = 0;
        labels.dt = 1;        
    end
    
    if ~isfield(labels,'tmax'), labels.tmax = labels.time(end); end
    if ~isfield(labels,'dt'), labels.dt = 1; end
    if ~isfield(labels,'stride'), labels.stride = 0; end
    if ~isfield(labels,'t0'), labels.t0 = 0; end
    
    if ~exist('commands','var')
        commands = '';
    end
    
    if isempty(xax), xax = 1:s(1); end;
    if isempty(yax), yax = 1:s(2); end;
       
    %% processing
  
    if stop == 1, 
        warning('Only one time step.');
        plotdata = double(squeeze(data)); % shiftdim screws up single snapshots
    else        
        plotdata = double(squeeze(shiftdim(data,index)));
    end    
    
    if labels.stride == 0
        datamax = nanmax(plotdata(:));
        datamin = nanmin(plotdata(:));
    else
        clim = caxis;
    end
    
    hfig = gcf;
    ckey = '';
    button = [];
    pflag = 0;
    spaceplay = 1; % if 1, space pauses. if 0, space plays
    
    %% parse options
    flag = [0 0 0 0 0 0];% defaults
    
    cmds = {'nocaxis','pcolor','imagesc','contour','pause','fancy_cmap'};
    for i = 1:length(cmds)
        loc = strfind(commands,cmds{i});
        if ~isempty(loc)
            flag(i) = 1;
            commands = [commands(1:loc-1) commands(loc+length(cmds{i}):end)];
        end
    end
    
    plotflag = sum([2 3 4] .* flag(2:4));
    if flag(5) && stop ~= 1, spaceplay = 0; fprintf('\n Hit a key to advance frame. \n\n'); end
    
    if flag(6) % Build colormap
        radius = 50;
        num = 40;
        theta = linspace(0, pi/2, num).';
        a = radius * cos(theta);
        b = radius * sin(theta);
        L = linspace(0, 100, num).';
        Lab = [L, a, b];
        fancy_map = applycform(Lab, makecform('lab2srgb'));
    end

    set(gcf,'Renderer','zbuffer'); % performance!
    i=0;
    while i<=stop-1

        % navigation
        if strcmp(ckey,'space') && isempty(button)
            spaceplay = ~spaceplay;
        end

        pflag = ~spaceplay;
        
        if pflag,
            [x,y,button] = ginput(1);
            figure(gcf);
            if button == 32, spaceplay = 1; end % resumes when paused
            if button == 27, break; end % exit when Esc is pressed.
        else
            pause(0.02);%(pausetime);
        end  
        
        ckey = get(gcf,'currentkey');% end
        
        % navigate : other keys move forward
        if strcmp(ckey,'leftarrow') | strcmp(ckey,'downarrow') | button == 28 | button == 31 | button == 3
            pflag = 1; spaceplay = 0;
            i = i-2;
        else if strcmp(ckey,'rightarrow') | strcmp(ckey,'uparrow') | button == 29 | button == 30 | button == 1
                pflag = 1; spaceplay = 0;
            end
        end
        
        if strcmp(ckey,'escape')
            break
        end
        
        i=i+1;
        if i < 1, i = 1; end
        
        %% Plot
        hold off; % just in case
%         switch plotflag
%             case 2
%                 pcolor(xax,yax,plotdata(:,:,i)');
%             case 3
%                 try
%                     imagescnan(yax,xax,plotdata(:,:,i)');
%                     set(gca,'yDir','normal');
%                 catch ME
%                     imagesc(xax,yax,plotdata(:,:,i)');
%                     set(gca,'yDir','normal')
%                 end
%             case 4
%                 set(gcf,'Renderer','painters');
%                 [C,h] = contour(xax,yax,plotdata(:,:,i)', 20,'k');
%                 format short
%                 clabel(C,h,'FontSize',9);
%             otherwise
%                 contourf(xax,yax,plotdata(:,:,i)', 40);
%         end
%         
%         % colorbar
%         if plotflag ~=4 
%             if ~flag(1)
%                 if labels.stride > 0
%                     caxis(clim);
%                 else
%                     if datamax ~= datamin, caxis([datamin datamax]); end
%                 end
%             end
%         end        
%         shading flat;
%         if flag(6), colormap(fancy_map); end
%         if plotflag ~= 4, colorbar;  end

        
        
        
        % labels
        if labels.revz, revz; end;
        if isempty(labels.time)
            title([labels.title ' t instant = ' num2str(labels.t0+i+(labels.dt-1)*(i-1)+100*labels.stride)]);
        else
            title([labels.title ' t = ' sprintf('%.2f/%.2f ', labels.time(i),labels.tmax) ...
                   labels.tunits ' (instant = ' num2str(labels.t0+i+(labels.dt-1)*(i-1)+100*labels.stride) ')']);
        end
        xlabel(labels.xax);
        ylabel(labels.yax);
        eval(commands); % execute custom commands
    end