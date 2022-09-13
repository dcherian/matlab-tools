% I want to fit y = y_0 tanh((x-x0)/X) + y_1
%  [y0, X, x0, y1, conf, fitobj] = tanh_fit(x, y, plot_flag, test)
function [y0, X, x0, y1, conf, fitobj] = tanh_fit(x, y, plot_flag, test)

    x = double(x); y = double(y);

    if ~exist('test', 'var'), test = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if test
        test_fit();
        return;
    end

    if size(x,1) ~= 1, x = x'; end
    if size(y,1) ~= 1, y = y'; end

    ym = mean(y);
    %warning('subtracting mean for better fit');
    y = y - ym;

    initGuess(1) = max(y(:));
    initGuess(2) = mean(x(:));
    initGuess(3) = mean(x(:));
    initGuess(4) = 0;

    opts = optimset('MaxFunEvals',1e7);
    [fit2,~,exitflag] = fminsearch(@(f) fiterror(f,x,y), ...
                                   initGuess,opts);

    y0 = fit2(1);
    X = fit2(2);
    x0 = fit2(3);
    y1 = fit2(4) + ym;

    tanhfn = fittype(@(y0,X,x0,y1,x) y0*tanh((x-x0)/X) + y1);
    fitobj = fit(x', y', tanhfn, 'StartPoint', fit2);
    conf = confint(fitobj);

    y0 = fitobj.y0;
    X = fitobj.X;
    if abs(conf(1,2)) > abs(conf(2,2))
        temp = conf(2,2);
        conf(2,2) = conf(1,2);
        conf(1,2) = temp;
    end
    x0 = fitobj.x0;
    y1 = fitobj.y1 + ym;
    conf(:,4) = conf(:,4) + ym;

    if ~exitflag
        fit2 = nan(4,1);
        conf(:) = NaN;
    end

    if plot_flag
        figure; hold on;
        plot(x,y,'k.', 'MarkerSize', 12);
        ci = predint(fitobj, x, 0.95, 'functional');
        plot(fitobj);
        plot(x, ci, 'r--');
        plot(fitobj,x',y', 'residuals');
        legend('off'); %'Location', 'NorthWest');
    end
end

function [E] = fiterror(f,x,y)
    y0 = f(1); X = f(2); x0 = f(3); y1 = f(4);

    E = sum((y - y0 .* tanh((x-x0)/X) - y1).^2);
end

function [] = test_fit()
    x = [0.02:0.05:4];
    X = 2;
    y0 = 2;
    x0 = 0;
    y1 = 2;
    y = y0 * tanh((x-x0)/X) + y1; % + rand(size(x)).*y0/4;

    [Y0, XX, X0, Y1] = tanh_fit(x,y,1);

    disp(['y0 = ' num2str(Y0) ' | Original = ' num2str(y0)]);
    disp(['X = ' num2str(XX) ' | Original = ' num2str(X)]);
    disp(['x0 = ' num2str(X0) ' | Original = ' num2str(x0)]);
    disp(['y1 = ' num2str(Y1) ' | Original = ' num2str(y1)]);
end