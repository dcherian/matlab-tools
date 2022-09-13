% Test coherence code

x = randn(1,10000);
n = randn(1,10000);

x = x/max(x);
n = n/max(n);

alpha = 6;
y = alpha*x+n;

[Cxy,F] = mscohere(x,y);
[fc,Cdc,pha] = dccoher(x,y,1,3,1);
Ca = alpha/sqrt(alpha^2 + var(n)/var(x));

fprintf('\n alpha = %.2f, Cal = %.4f, MATLAB = %.4f, DC = %.4f \n', alpha,Ca,mean(Cxy), mean(abs(Cdc).^2));