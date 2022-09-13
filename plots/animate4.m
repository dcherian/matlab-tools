% Makes animation (using contourf). Assumes last dimension is to be looped by default. 
% Else specify index. Allows browsing.
%       [] = animate(data,labels,commands,index,pausetime)
%           xax,yax - CELLS with 4 components each
%           data - CELL - data to be animated and axes - script squeezes data
%           labels - CELL - Structures with following fields for labeling plot (optional - can be [])
%               > title
%               > xax
%               > yax
%               > revz - reverse yDir?
%               > time - where to pull time info for titles
%               > tunits - units for time vector labels.time
%           commands - custom commands to execute after plot or one of below. (in string, optional)
%                    - separate commands using ;
%                    - APPPLIES TO ALL PLOTS
%               > nocaxis - non-constant colorbar
%           index - dimension to loop through (optional)
%           pausetime - pause(pausetime) (optional)
%
% USAGE:
%       animate(data)
%       animate(xax,yax,data)
%       animate([],[],data,...
%       animate(xax,yax,data,commands)
%
% BROWSE:
%       - first space pauses, second space resumes, remaining spaces *play*
%       - arrowkeys *pause* and navigate always
%

function [] = animate4(xax,yax,data,labels,commands,index,pausetime)

    % figure out input
    narg = nargin;
    
%     switch narg
%         case 1,
%             data = squeeze(xax);
%             s    = size(data);
%             xax  = 1:s(1); 
%             yax  = 1:s(2);
%             
%         case 4,
%             if strcmp(class(labels),'char')
%                commands = labels;
%                labels = [];
%             end
%     end
    
    data = squeeze(data);
    s = size(data);
       
    if narg <= 4 || ~exist('index','var')
        stop = size(data,3);
        index = length(s); % set last dimension to loop by default
    else
        stop = s(index);
    end
    
    if narg ~= 6 || ~exist('pausetime','var')
        pausetime = 0.2;
    end
    
    if ~exist('labels','var') || isempty(labels)
        labels.title = '';
        labels.xax   = '';
        labels.yax   = '';
        labels.revz  = 0;
        labels.time  = [];
    end
    
    if ~exist('commands','var')
        commands = '';
    end
    
    if isempty(xax), xax = 1:s(1); end;
    if isempty(yax), yax = 1:s(2); end;
    
    %% processing
    
    if stop == 1, warning('Only one time step.'); end
    
    plotdata{1} = double(squeeze(shiftdim(data{1},index)));
    plotdata{2} = double(squeeze(shiftdim(data{2},index)));
    plotdata{3} = double(squeeze(shiftdim(data{3},index)));
    plotdata{4} = double(squeeze(shiftdim(data{4},index)));
    
    datamax{1} = nanmax(plotdata{1}(:));
    datamin{1} = nanmin(plotdata{1}(:));
    datamax{2} = nanmax(plotdata{2}(:));
    datamin{2} = nanmin(plotdata{2}(:));
    datamax{3} = nanmax(plotdata{3}(:));
    datamin{3} = nanmin(plotdata{3}(:));
    datamax{4} = nanmax(plotdata{4}(:));
    datamin{4} = nanmin(plotdata{4}(:));
    
    hfig = figure;
    ckey = '';
    button = [];
    pflag = 0;
    spaceplay = 1; % if 1, space pauses. if 0, space plays
    caxisflag = 1; % constant color bar
    
    loc = strfind(commands,'nocaxis');
    if ~isempty(loc)
        caxisflag = 0;
        commands = [commands(1:loc-1) commands(loc+length('nocaxis'):end)];
    end
    
    i=0;
    while i<=stop-1

        if strcmp(ckey,'space') & isempty(button)
            spaceplay = ~spaceplay;
        end

        pflag = ~spaceplay;
        
        if pflag,
            [x,y,button] = ginput(1);
            figure(gcf);
            if button == 32, spaceplay = 1; end % resumes when paused
        else
            pause(pausetime);
        end  
        
        ckey = get(gcf,'currentkey'); 
        
        % navigate : other keys move forward
        if strcmp(ckey,'leftarrow') | strcmp(ckey,'downarrow') | button == 28 | button == 31 
            pflag = 1; spaceplay = 0;
            i = i-2;
            if i < 1, i = 1; end
        else if strcmp(ckey,'rightarrow') | strcmp(ckey,'uparrow') | button == 29 | button == 30 
                pflag = 1; spaceplay = 0;
            end
        end
        
        i=i+1;
        
        %% Plot
        
        for i=1:4
            plotvar = plotdata{i};
            contourf(xax{i},yax{i},plotvar(:,:,i)', 40);
            if labels{i}.revz, revz; end;
            shading flat
            if isempty(labels{i}.time)
                title([labels.title{i} ' t instant = ' num2str(i)]);
            else
                title([labels.title{i} ' t = ' sprintf('%.2f/%.2f ', labels.time(i),labels.time(end)) labels.tunits{i}]);
            end
            xlabel(labels.xax{i});
            ylabel(labels.yax{i});
            if caxisflag
                if datamax{i} ~=datamin{i}, caxis([datamin{i} datamax{i}]); end
            end
            colorbar;        
            eval(commands); % execute custom commands
        end
    end
    
    
        % debug
        %fprintf('\n i = %2d pflag = %d button = %3d ckey = %5s', i, pflag, button, ckey);
        
%         if spaceplay
%             pflag = 0;
%         else
%             pflag = 1;
%         end