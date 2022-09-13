function [] = TestBootstrap

% generate noisy time Series
    x = [0:1/8000:1]';
    y = x + 0.1*rand(size(x));

    plot(x,y,'*');

    [b, bint] = regress(y,x);
    bootr = @(x,y)regress(y,x);

     se = bootci(100, {bootr, x, y}, 'type', 'student');