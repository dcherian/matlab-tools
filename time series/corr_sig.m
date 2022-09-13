function y=significance(N,sig_level)

% routine compute correlation coefficient corresponding to sig_level%
% significance level for N degrees of freedom

% for example, to find the correlation coefficient which has a 95%
% chance of being different from zero---that is in only 5% of cases
% computed with random noise would you get a correlation coefficient
% greater than y.
%   y=significance(200,.95);

y=erfinv(sig_level)*sqrt(2/N);

% to test this out with white noise:
%xx=randn(100,100);
%yy=randn(100,1);
%for i=1:100
% rr=corrcoef(yy,xx(:,i));
% zz(i)=rr(2,1);
%end
%zz2=sort(abs(zz));

% compare:  zz2(50) with erfinv(.5)*sqrt(2/100)

