function [] = probdens(x)

edges = linspace(min(x(:)),max(x(:)),20);
N = histc(x(:),edges);
bar(edges,N./numel(x),'histc')
xlim([min(x(:)) max(x(:))])