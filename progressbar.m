% To use : 
%                 cpb = progressbar()
%        in loop: progressbarupdate(cpb,VALUE,TEXT); - VALUE in %
%     after loop: cpb.stop();

function [cpb] = progressbar()

cpb = ConsoleProgressBar();

% Set progress bar parameters
cpb.setLeftMargin(1); % progress bar left margin
cpb.setTopMargin(1); % rows margin

cpb.setElapsedTimeVisible(1)
%cpb.setRemainedTimeVisible(1)
cpb.setElapsedTimePosition('left')
%cpb.setRemainedTimePosition('right')

cpb.setLength(40); % progress bar length: [.....]
cpb.setMinimum(0); % minimum value of progress range [min max]
cpb.setMaximum(100); % maximum value of progress range [min max]

cpb.start();