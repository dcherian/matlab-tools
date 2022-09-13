function [q] = gcm_pv(file,timesteps)

    file = find_file(file);   
    
    if  ~exist('timesteps','var') || isempty(timesteps)
        timesteps = [1 Inf];
        dt = 1;
    end
    
    switch(numel(timesteps))
        case 1
            timesteps(2) = timesteps(1);
            dt = 1;
        case 2
            dt = 1;
        case 3
            dt = timesteps(2);
            timesteps = [timesteps(1) timesteps(2)];
    end
    
    f = 1E-4;
    Xp1 = ncread(file,'Xp1');
    Yp1 = ncread(file,'Yp1');
    X   = ncread(file,'X');
    Y   = ncread(file,'Y');
    Z   = ncread(file,'Z');
    T   = ncread(file,'T');
    
    if isinf(timesteps(2))
        T = T(timesteps(1):dt:end);
    else
        T = T(timesteps(1):dt:timesteps(2));
    end
    
    u    = ncread(file,'U',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1) + 1)]);
    v    = ncread(file,'V',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1)+ 1)]);
    temp = ncread(file,'Temp',[1 1 1 timesteps(1)],[Inf Inf Inf (timesteps(2)-timesteps(1)+ 1)]); % theta
    
    lambda = temp;
    
    s = size(lambda);
    if length(s) == 3
        s(4) = 1;
    end    
    
    [xu, yu, zu, tu] = ndgrid(Xp1,Y,Z,T);
    [xv, yv, zv, tv] = ndgrid(X,Yp1,Z,T);
    [xl, yl, zl, tl] = ndgrid(X,Y,Z,T);
    
    fmat = repmat(f,size(lambda));
   
    %q = pv(Xp1,Y,Z,u, X,Yp1,Z,v, X,Y,Z,temp,f);

    vx    = diff(v,1,1)./diff(xv,1,1);%repmat(diff(X,1,1),[1 s(2)+1 s(3) s(4)]);
    vy    = diff(v,1,2)./diff(yv,1,2);%repmat(diff(Yp1',1,2),[s(1) 1 s(3) s(4)]);
    vz    = diff(v,1,3)./diff(zv,1,3);%repmat(diff(Z,1,1) ,[s(1) s(2)  1   s(4)]);

    ux    = diff(u,1,1)./diff(xu,1,1);%repmat(diff(xu',1,1),[1 1 s(3) s(4)]);
    uy    = diff(u,1,2)./diff(yu,1,2);%repmat(diff(yu',1,2),[1 1 s(3) s(4)]);
    uz    = diff(u,1,3)./diff(zu,1,3);%repmat(diff(zu,1,1) ,[1 1  1  s(4)]);

    lx    = diff(lambda,1,1)./diff(xl,1,1);%repmat(diff(xl',1,1),[1 1 s(3) s(4)]);
    ly    = diff(lambda,1,2)./diff(yl,1,2);%repmat(diff(yl',1,2),[1 1 s(3) s(4)]);
    lz    = -diff(lambda,1,3)./diff(zl,1,3);%repmat(diff(zl,1,1) ,[1 1  1   s(4)]);
    
    index = find(s(1:3) == 1); % find singleton dimensions and correct derivative (exclude time)
    switch index
        case 1
            %ux = zeros(size(xu)-[1 0 0 0]);
            vx = zeros(size(xv));%-[1 0 0 0]);
            lx = zeros(size(xl));%-[1 0 0 0]);
        case 2
            uy = zeros(size(xu));%-[0 1 0 0]);
            %vx = zeros(size(xv)-[1 0 0 0]);
            ly = zeros(size(xl));%-[0 1 0 0]);
        case 3
            uz = zeros(size(xu));%-[0 0 1 0]);
            vz = zeros(size(xv));%-[0 0 1 0]);
            lz = zeros(size(xl));%-[0 0 1 0]);
    end
    
    coeff1 = (fmat(1:end-1,:,:,:) + vx(:,1:end-1,:,:) - uy(2:end-1,:,:,:));
    pv1    = coeff1(:,:,1:end-1,:).*lz(1:end-1,:,:,:);
    pv2    = (-1)*vz(:,1:end-1,:,:).*ly(:,:,1:end-1,:);
    pv3    = uz(2:end-1,:,:,:).*lx(:,:,1:end-1,:); 
    q      = (pv1 + pv2(1:end-1,:,:,:) + pv3)./lambda(1:end-1,:,1:end-1,:);

    %% for debugging purposes
    
    %     
%     Vx    = diff(v,1,1)./diff(xv,1,1);%repmat(diff(X,1,1),[1 s(2)+1 s(3) s(4)]);
%     Vy    = diff(v,1,2)./diff(yv,1,2);%repmat(diff(Yp1',1,2),[s(1) 1 s(3) s(4)]);
%     Vz    = diff(v,1,3)./diff(zv,1,3);%repmat(diff(Z,1,1) ,[s(1) s(2)  1   s(4)]);
% 
%     Ux    = diff(u,1,1)./diff(xu,1,1);%repmat(diff(xu',1,1),[1 1 s(3) s(4)]);
%     Uy    = diff(u,1,2)./diff(yu,1,2);%repmat(diff(yu',1,2),[1 1 s(3) s(4)]);
%     Uz    = diff(u,1,3)./diff(zu,1,3);%repmat(diff(zu,1,1) ,[1 1  1  s(4)]);
% 
%     Lx    = diff(lambda,1,1)./diff(xl,1,1);%repmat(diff(xl',1,1),[1 1 s(3) s(4)]);
%     Ly    = diff(lambda,1,2)./diff(yl,1,2);%repmat(diff(yl',1,2),[1 1 s(3) s(4)]);
%     Lz    = -diff(lambda,1,3)./diff(zl,1,3);%repmat(diff(zl,1,1) ,[1 1  1   s(4)]);
