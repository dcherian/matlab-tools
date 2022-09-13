% I want to fit y = y_0 exp(-((x-x0)/X)^2) + y1
%      [y0, X, x0, y1, exitflag, conf] = gauss_fit(x, y, plot_flag, test)
function [y0, X, x0, y1, exitflag, conf] = gauss_fit(x, y, plot_flag, test)

    x = double(x); y = double(y);

    if ~exist('test', 'var'), test = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if nargin == 0
        disp('Running test')
        test = 1;
    end

    if test
        test_fit();
        return;
    end

    if size(x,1) ~= 1, x = x'; end
    if size(y,1) ~= 1, y = y'; end

    initGuess(1) = max(y);
    [initGuess(2), index] = max(x(:));
    initGuess(3) = x(index);
    initGuess(4) = min(abs(y(:)));
    opts = optimset('MaxFunEvals',1e11);

    %this doesn't seem to work
    %modelfun = @(b,t)b(1)*exp(-(t-b(2))/b(3)) + b(4);
    %fitmodel = fitnlm(x,y, modelfun, initGuess);

    % do the fit using fminsearch.
    % use solution to curve-fit so that MATLAB gives me confidence intervals :)
    % launda isshmart hai!
    gaussfit = fittype(@(y0,X,x0,y1,x)y0*exp( -((x-x0)/X).^2 ) + y1);
    [fit2,~,exitflag] = fminsearch(@(fit) fiterror(fit,x,y), ...
                                   initGuess,opts);
    fitobj = fit(x',y',gaussfit, 'StartPoint', fit2);
    conf = confint(fitobj);

    y0 = fitobj.y0;
    X = abs(fitobj.X); % sometimes returns -ve for some reason
    conf(:,2) = abs(conf(:,2));
    if conf(1,2) > conf(2,2)
        temp = conf(2,2);
        conf(2,2) = conf(1,2);
        conf(1,2) = temp;
    end
    x0 = fitobj.x0;
    y1 = fitobj.y1;

    if plot_flag
        figure;
        subplot(211);
        hold on;
        plot(x,y,'k.', 'MarkerSize', 12);
        plot(fitobj, 'predobs');
        plot(fitobj,x',y', 'residuals')
    end
end

function [E] = fiterror(fit,x,y)
% x = (T0,H,a)
    y0 = fit(1); X = fit(2); x0 = fit(3); y1 = fit(4);

    E = sum((y - y0 .* exp(-((x-x0)/X).^2) - y1).^2);
end

function [] = test_fit()
    x = [-24:0.05:24];
    X = 10;
    y0 = 2;
    x0 = 0.5;
    y1 = 0.2;
    y = y0 * exp(-((x-x0)/X).^2) + y1 + y0/50 * rand(size(x));

    [yy,xx,xx0, y10] = gauss_fit(x,y,1);

    disp(['y0 = ' num2str(yy) ' | Original = ' num2str(y0)]);
    disp(['X = ' num2str(xx) ' | Original = ' num2str(X)]);
    disp(['x0 = ' num2str(xx0) ' | Original = ' num2str(x0)]);
    disp(['y1 = ' num2str(y10) ' | Original = ' num2str(y1)]);
end
