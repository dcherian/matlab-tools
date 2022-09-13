% Does Ordinary least squares. Calculates confidence bounds after
% estimating degrees of freedom.
%    [coeff,conf,dof,err] = dcregress(x, y, dof, test_flag, plot_flag, dof_flag)
% coeff = [intercept slope]
% conf = 95% confidence intervals on coeff
% dof =  degrees of freedom used for confidence interval
% err = standard error on coeff.
% dof_flag = 1 => if dof not provided estimate using calcdof
%          = 0 => all input points are independent
% err_flag = 1 => calculate errors else skip it.
function [coeff,conf,dof,err] = dcregress(x, y, dof, test_flag, plot_flag, ...
                                          dof_flag, err_flag)

    if ~exist('dof', 'var') | isempty(dof), dof = NaN; end
    if ~exist('test_flag', 'var'), test_flag = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end
    if ~exist('err_flag', 'var'), err_flag = 1; end
    if ~exist('dof_flag', 'var'), dof_flag = 1; end

    if test_flag
        test_dcregress;
        return;
    end

    if size(x,2) == 1, x = x'; end
    if size(y,2) == 1, y = y'; end

    xnan = x; ynan = y;
    nanmask = isnan(x) | isnan(y);
    x = xnan(~nanmask);
    y = ynan(~nanmask);

    E = [ones(size(x')) x'];

    %%%%%%%%%%% See Wunsch(1996) pg. 116
    % P matrix
    coeff = E\y';
    intercept = coeff(1);
    slope = coeff(2);

    if ~err_flag & ~plot_flag, return; end
    true = y; est = intercept + slope .* x;
    res = true-est;
    % (E' * E) ^-1
    %ETEI = inv(E'*E);
    % from http://blogs.mathworks.com/loren/2007/05/16/purpose-of-inv/
    [Q,R] = qr(E,0);
    S = inv(R);
    ETEI = S*S';
    % assuming noise vector (res) is white
    P = ETEI * E' * var(res) * E * ETEI;
    err = sqrt(diag(P));  % standard error

    % Rnn = res'*res;
    % durbin test
    % rhohat = sum(res(2:end).*res(1:end-1)) ...
    %    ./ sum(res.^2)

    %assert((seb1-err(2))./err(2) < 1e-2);

    if isnan(dof)
        if dof_flag
            % calculate degrees of freedom taking into account gaps
            dof(1) = calcdof(xnan);
            dof(2) = calcdof(ynan);
            dof = min(dof-2);
        else
            % assume all points are statistically independent
            dof = length(x)-2;
        end
    end
    % standard error of slope assuming uncorrelated Gaussian white noise.
    % n = dof; % length(x);
    % r = corrcoef(x,y); r = r(1,2);
    % syx = ((n-1)/(n-2) * var(y) * (1-r^2))^(1/2);
    % seb1 = syx/sqrt(var(x)*(n-1)); % should equal err(2)

    [lo,up] = conft(0.05, dof);
    conf(1,1) = up * err(1);
    conf(2,1) = up * err(2);

    %%%%%%%%%%% use MATLAB regress
    %[b, bint, r, rint, stats] = regress(y', E);

    % slope = b(2);
    % err = abs(bint(2) - b(2));

    % plot fit
    % figure; hold all;
    % plot(tvec/86400, true, '*');
    % plot(tvec/86400, est); plot(tvec/86400, res); liney(0);

    %[c,lags] = xcorr(fluxvec - mean(fluxvec(:)), 'coef');
    %plot(lags, c); linex(0); liney(0);

    if plot_flag
        figure;
        subplot(2,3,[1 2]);
        hold on;
        plot_fit(x, y, coeff(2), conf(2), coeff(1), conf(1), 'k', 1);
        title(['Slope = ' num2str(coeff(2)) ' | ' ...
              'Intercept = ' num2str(coeff(1))])
        %rr = mf_wtls(x,y,0,0.4*std(y));
        %plot_fit(x, y, rr(1), rr(2), rr(3), rr(4), 'r', 0);
        xlabel('x'); ylabel('y');
        beautify;

        hax(1) = subplot(234);
        histogram(res, 40, 'Normalization', 'PDF');
        hold on;
        xlim([-1 1]*max(abs(xlim)));
        %FitAndPlotDist(res, 'normal', hax(1));
        xplt = linspace(min(xlim), max(xlim), 100);
        hnorm = plot(xplt, normpdf(xplt, 0, std(res)), 'b-');
        legend(hnorm, 'Normal \mu=0, \sigma=std(res)')
        ylabel('residuals PDF')
        beautify;

        subplot(235);
        [c, lags] = xcorr(res, 'coef');
        plot(lags, c);
        hold on;
        xlim([-100 100])
        liney(0); linex(0);
        ylabel({'Autocorrelation'; 'of residuals'})
        beautify;

        subplot(233)
        plot(x, res, '*');
        xlabel('x'); ylabel('residuals'); liney(0);
        beautify;

        subplot(236)
        cla
        n = 1;
        plot(res(1:end-n), res(n+1:end), '*');
        xlabel('residuals(t)');
        ylabel(['residuals(t+' num2str(n) ')']);
        liney(0); linex(0);
        beautify;
    end
end

function [] = plot_fit(x, y, slp, slperr, int, interr, color, plotxy)

    if plotxy
        plot(x,y, '.', 'MarkerSize', 16, 'Color', [1 1 1]*0.5);
    end
    limx = xlim;
    xvec = linspace(limx(1), limx(2), 1000);
    yvec = slp * xvec + int;
    yvechi = (slp + slperr)*xvec + int + interr;
    yveclo = (slp - slperr)*xvec + int - interr;

    plot(xlim,xlim.*slp+int,'k--','LineWidth',2)
    plot([0 max(xlim)],[0 max(xlim)].*(slp+slperr)+int+interr,'k--', ...
         'Color', color, 'LineWidth',0.75)
    plot([0 max(xlim)],[0 max(xlim)].*(slp-slperr)+int-interr,'k--', ...
         'Color', color, 'LineWidth',0.75)
    plot([min(xlim) 0],[min(xlim) 0].*(slp+slperr)+int-interr,'k--', ...
         'Color', color, 'LineWidth',0.75)
    plot([min(xlim) 0],[min(xlim) 0].*(slp-slperr)+int+interr,'k--', ...
         'Color', color, 'LineWidth',0.75)
end

function [] = test_dcregress()

    b = 0.1
    a = 0.5
    x = linspace(0,1,8000);
    y = a*x + b + a/4 * randn(size(x));

    [coef,conf] = dcregress(x,y, [], 0, 1);
    [coef-conf coef+conf]
    % [b, bint, r, rint, stats] = regress(y', [ones(size(x')) x']);
end
