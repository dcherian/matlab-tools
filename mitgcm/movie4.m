function [] = movie4(file,timesteps)

    file = find_file(file);
    
    if ~exist('timesteps','var') || isempty(timesteps)
        timesteps = [1 Inf];
        dt = 1;
    else    
        if length(timesteps) == 3
            dt = timesteps(2);
            timesteps(2) = timesteps(3);
            timesteps(3) = NaN;
        else
            if length(timesteps) == 2
                dt = 1;
            end    
        end
    end
    
    Xp1 = ncread(file,'Xp1');
    Yp1 = ncread(file,'Yp1');
    X   = ncread(file,'X');
    Y   = ncread(file,'Y');
    Z   = ncread(file,'Z');
    T   = ncread(file,'T',timesteps(1),timesteps(2),dt);
    
    stride = [1 1 1 dt];
    
    u    = ncread(file,'U',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1) + 1)], stride);
    v    = ncread(file,'V',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1)+ 1)], stride);
    w    = ncread(file,'W',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1)+ 1)], stride);
    temp = ncread(file,'Temp',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1)+ 1)], stride); 
    
    %ulim = [min(u(:)) max(u(:))];
    
    s = size(u);
    s1 = size(v);
    
    if repnan(u(1,1,:,:),0)   == zeros([1 1 s(3) s(4)]), u(1,1,:,:) = NaN; end;
    if repnan(u(end,1,:,:),0) == zeros([1 1 s(3) s(4)]), u(end,1,:,:) = NaN; end;
    if repnan(v(1,1,:,:),0)   == zeros([1 1 s(3) s(4)]), v(1,1,:,:) = NaN; end;
    if repnan(v(end,1,:,:),0) == zeros([1 1 s(3) s(4)]), v(end,1,:,:) = NaN; end;
    if repnan(temp(1,1,:,:),0)   == zeros([1 1 s(3) s(4)]), temp(1,1,:,:) = NaN; end;
    if repnan(temp(end,1,:,:),0) == zeros([1 1 s(3) s(4)]), temp(end,1,:,:) = NaN; end;
    if repnan(w(1,1,:,:),0)   == zeros([1 1 s(3) s(4)]), w(1,1,:,:) = NaN; end;
    if repnan(w(end,1,:,:),0) == zeros([1 1 s(3) s(4)]), w(end,1,:,:) = NaN; end;
    
    stop=size(u,4);
    
    hfig = figure;
    ckey = '';
    button = [];
    pflag = 0;
    spaceplay = 1; % if 1, space pauses. if 0, space plays
    %caxisflag = 1; % constant color bar
    
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
            if button == 27, break; end % exit when Esc is pressed.
        %else
            %pause(pausetime);
        end  
        
        ckey = get(gcf,'currentkey'); 
        
        % navigate : other keys move forward
        if strcmp(ckey,'leftarrow') | strcmp(ckey,'downarrow') | button == 28 | button == 31 
            pflag = 1; spaceplay = 0;
            i = i-2;
        else if strcmp(ckey,'rightarrow') | strcmp(ckey,'uparrow') | button == 29 | button == 30 
                pflag = 1; spaceplay = 0;
            end
        end
        
        if strcmp(ckey,'escape')
            break
        end
        
        i=i+1;
        if i < 1, i = 1; end
        
        %% Plot
        %plotflag
        
        subplot(221)
    	contourf(X,Z,squeeze(v(:,1,:,i+1))',40);
        xlabel('X')
        ylabel('Z')
        title('V');
        shading flat; colorbar
        
        subplot(222)
    	contourf(Xp1,Z,squeeze(u(:,1,:,i+1))',40);
        xlabel('X')
        ylabel('Z')
        title('U');
        shading flat; colorbar
        
        subplot(223)
    	contourf(X,Z,squeeze(temp(:,1,:,i+1))',40);
        xlabel('X')
        ylabel('Z')
        title('Temp');
        shading flat; colorbar
        
        subplot(224)
    	contourf(X,Z,squeeze(w(:,1,:,i+1))',40);
        xlabel('X')
        ylabel('Z')
        title('W');
        shading flat; colorbar
        
        suplabel(sprintf('t = %.2f/%.2f days', (double(T(i))/3600/24), double(T(end))/3600/24), 't');
        
   %     colorbar;        
   %     eval(commands); % execute custom commands
    end
    
    
    
    