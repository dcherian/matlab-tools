% I want to fit y = y0*(1 - erf(x/X))
function [y0, X] = erf_fit(x, y, plot_flag, test)

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

    initGuess(1) = max(y);
    initGuess(2) = x(ceil(length(x)/2));
    initGuess(3) = 0;
    opts = optimset('MaxFunEvals',1e7);
    [fit2,~,exitflag] = fminsearch(@(fit) fiterror(fit,x,y), ...
                                   initGuess,opts);

    y0 = fit2(1);
    X = (fit2(2)); % sometimes returns -ve for some reason

    if plot_flag
        figure;
        plot(x,y,'k*'); hold all
        plot(x, y0*(1-erf(x/X)));
        legend('Data','Fit')
    end
end

function [E] = fiterror(fit,x,y)
% x = (T0,H,a)
    y0 = fit(1); X = fit(2);

    E = sum((y - y0 .* (1-erf(x/X))).^2);
end

function [] = test_fit()
    x = [0:0.05:3];
    X = 1;
    y0 = -2;
    y = y0 * (1-erf(x/X)); % + y0/100 * rand(size(x));

    [yy,xx] = erf_fit(x,abs(y),1);

    disp(['y0 = ' num2str(yy) ' | Original = ' num2str(y0)]);
    disp(['X = ' num2str(xx) ' | Original = ' num2str(X)]);
end
