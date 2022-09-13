function [Vmode, Hmode, c] = vertmode(N2, Z, n, make_plot, check_norm)

% function [Vmode, Hmode, c] = vertmode(N2, Z, n, make_plot, check_norm)
% Uses centered first differences discretization of Chelton et al. (1998)
% Takes input N2 -> Buoyancy frequency squared (@ mid pts of Z)
%              Z -> Vertical co-ordinate
%              n -> No. of modes to isolate
%                   Max. number of nodes is length(Z)-1
%              make_plot -> plot mode shapes
%              check_norm -> check normalization
% Returns
%             Vmode -> Vertical structure for vertical velocity
%             Hmode -> Vertical structure for horizontal velocities, pressure
%              c(i) -> Gravity Wave speed of i-th mode
%

    if ~exist('make_plot', 'var'), make_plot = 1; end
    if ~exist('check_norm', 'var'), check_norm = 0; end

    lz = length(Z);
    if size(Z,1) == 1, Z = Z'; end
    if size(N2,1) == 1, N2 = N2'; end

    assert(length(N2) == length(Z)-1);
    if any(N2 < 0)
       error('N² < 0, won''t work!');
    end

    if all(Z < 0)
        warning('Provided negative Z - converting to positive');
        Z = -1 * Z;
    end

    if any(Z < 0)
        error(['Some elements in Z are above the surface!']);
    end

    flipflag = 0;
    if all(diff(Z) < 0),
        warning('Flipping inputs so that diff(Z) > 0 always');
        N2 = flipud(N2);
        Z = flipud(Z);
        flipflag = 1;
    end

    if n > (lz-1)
        fprintf('\n n too big. Reducing to (length(Z) - 1)');
        n = lz - 1;
    end

    % Following Chelton et al. (1998)
    % Q'AW = (lambda)W (eqn. C.5)
    % centered first difference discretization
    % of the eigenvalue problem

    D(1) = Z(2)/2;
    D(2:lz-1) = (Z(3:lz)-Z(1:lz-2))/2;
    D(lz) = (Z(lz) - Z(lz-1))/2;
    Zmid = (Z(1:lz-1) + Z(2:lz))/2;

    Q1 = zeros(lz-1,1);

    for k = 1:lz-1
        Q1(k) = 2/(D(k)*D(k+1)*(D(k) + D(k+1)));
    end

    Q = diag(Q1./N2);
    A = -1*(diag(D(1:lz-1)) + diag(D(2:lz))) + diag(D(1:lz-2),1) + diag(D(3:lz),-1);

    [G,e] = eig(Q*A);

    c = sqrt(-1./diag(e));
    if ~isreal(c), error('Imaginary c!'); end
    [~,ind] = sort(c); % ascending order

    F=zeros(lz-2,lz-1);
    for i=1:lz-1
        F(:,i) = (c(i)^2/9.81)*(diff(G(:,i))./diff(Zmid));
    end

    Hmode = fliplr(F(:,ind(lz-n:lz-1)));
    Vmode = fliplr(G(:,ind(lz-n:lz-1)));

    % Fill in Hmode at lowest depth
    Hmode(lz-1,:) = Hmode(lz-2,:);
    c = flipud(c(ind(lz-n:lz-1)));

    % Normalize by max. amplitude
    %Hmode = Hmode./repmat(max(abs(Hmode)),length(Hmode),1);
    %Vmode = Vmode./repmat(max(abs(Vmode)),length(Vmode),1);

    % Normalize by energy (following Wunsch(1999)
    % why am I doing this using Hmode norm?
    norm = sqrt(sum(avg1(Vmode).^2.*repmat(diff(Zmid),1,n)));
    Vnorm = Vmode./repmat(norm,lz-1,1);
    norm = sqrt(sum(avg1(Hmode).^2.*repmat(diff(Zmid),1,n)));
    Hnorm = Hmode./repmat(norm,lz-1,1);

    if check_norm
        % check normalization
        hchk = sum(avg1(Hnorm).^2.*repmat(diff(Zmid),1,n));
        vchk = sum(avg1(Vnorm.*repmat(N2,1,n)).^2.*repmat(diff(Zmid),1,n));

        fprintf('\n u mode normalization: \n');
        disp(hchk);
        fprintf('\n w mode normalization: \n');
        disp(vchk);
    end

    Hmode = Hnorm;
    Vmode = Vnorm;

    if flipflag
        Hmode = flipud(Hmode);
        Vmode = flipud(Vmode);
    end

    % Plot first n modes
    if make_plot
        figure;
        ax(1) = subplot(131);
        plot(N2,Zmid);
        title('N^2');
        ylabel('Z (m)');
        if Z > 0, set(gca,'ydir','reverse'); end

        ax(2) = subplot(132);
        plot(Vmode,Zmid);
        hold on;
        linex(0);
        ylabel('Z (m)');
        title('Vertical Structure of Vertical Velocity');
        legend(num2str((1:n)'));
        if Z > 0, set(gca,'ydir','reverse'); end
        beautify;

        ax(3) = subplot(133);
        plot(Hmode,Zmid);
        hold on;
        linex(0);
        plot(xlim,[4000 4000], 'k--');
        ylabel('Z (m)');
        title('Vertical Structure of Horizontal Velocity');
        legend(num2str([1:n]'));
        if Z > 0, set(gca,'ydir','reverse'); end
        linkaxes(ax,'y');

        beautify;
    end
