% compute baroclinicity - better defined as metric of vertical
% uniformity of a profile.
%     bc = (KE_total - KE_depthavg)/KE_tot 
%          where KE =  ∫_z 0.5 * profile.^2 dz

function [bc] = baroclinicity(zvec, profile)

    if ~exist('zvec', 'var') && ~exist('profile', 'var')
        % test mode
        test_flag = 1;
        
        zvec = [0:0.005:1];
        profile = (sin(pi*zvec))/2;

        % analytic
        syms z
        prof = (sin(pi*z))/2
        mean = int(prof,[0 1]);
        tot = 0.5 * int(prof.^2, [0 1]);
        da = 0.5 * mean.^2 * 1;
        result = (tot - da)/tot;
    else
        test_flag = 0;
    end

    pmean = nanmean(profile);

    ketot = abs(0.5 * trapz(zvec, profile.^2));
    keda = abs(0.5 * trapz(zvec, pmean^2 .* ones(size(profile))));

    bc = (ketot - keda)/ketot;

    if test_flag
        figure;
        plot(profile, zvec);
        title(['bc = ' num2str(bc) ' v/s analytic = ' num2str(double(result))]);
    end