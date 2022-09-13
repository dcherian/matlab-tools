% -------------------------------------------------------------------------
%                       RANDOM
% ------------------------------------------------------------------------- 
% avg1                      - calculates average of all sucessive 2 point intervals
% beautify_session          - modifies various figure/axes/line/text object
%                             properties for the session.
% clr                       - clears all variables and closes all figures
%
% sdisp                     - Pretty print simplified symbolic expression
%
% animate                   - Shows movie of input data.s
% 
% check_gap                 - 
% fill_data                 - replaces NaN's with linearly interpolated values
% fill_gap                  -
% find_gap                  -
%
% simple_ls                 - simple least squares by right division
% taper_ls                  - smoothed & tapered least squares
% svd_fit                   - svd solution of system
%
% addnan                    - replaces values > input with NaN.
% fillnan                   - replaces values == input with NaN. More reasonable 
% repnan                    - Replaces NaN's with specified value
%
% gen_ser                   - generates a series 1:size(variable,dimension)
% nan_detrend               - Removes columnwise mean from each column of var.
% revz                      - reverses ydir
% stat                      - returns variable statistics. Accepts string input.
% txp                       - calculates 10^(). Accepts string input
%
% -------------------------------------------------------------------------
%                       OCEANOGRAPHY
% -------------------------------------------------------------------------
% bfrq                      - calculates buoyancy frequency
% vertmode                  - calculates vertical modes
%
% roms_movie                - displays movie of variable in ROMS output file
% roms_diffmovie            - displays movie of variable - variable at
%                             first "time" instant.
% roms_pv                   - calculates Ertel PV from ROMS output file
%
% -------------------------------------------------------------------------
%                       TIME SERIES STUFF
% -------------------------------------------------------------------------
% aliasfreq                 - returns apparent frequency given a known frequency
% confchi2                  - Confidence interval for chi-squared variate
% conft                     - Confidence interval for t-distributed variate
% dcdetide                  - removes the main 4 tidal constituents
% dcconv                    - calculates discrete convolution (1D only). replicates MATLAB's conv. use that instead.
% dccoher                   - Calculates coherence amplitude and phase.
% dcfft                     - calculates fft and returns power, frequency, FFT coeff's
% dchist2                   - plots 2d histogram (BETTER THAN BELOW)
% dchist2d                  - calls hist2d (in misc) with appropriate arguments and plots result.
% dcpgram                   - plots periodogram of input. does fft inside and marks top 'x' peaks
% dcpsd                     - computes and plots power spectrum density
% chkparsvl                 - takes fft coeffs and data and checks whether parseval's  theorem is satisfied. 
%                           - dcfft calls this by default.
%
% -------------------------------------------------------------------------
%                           PLOTS
% -------------------------------------------------------------------------
% beautify                  - call after making plot to make kickass
% datex                     - calls datetick on x-axis
% disp_plot                 - plots time series with depth displaced by value of depth
% loglinex(x,factor)        - Plots and labels vertical line at a given x and a factor 
%                             to scale the labelled 1./x value by.%
% fix_subplot2x2            - Tries to reduce extra whitespace in 2x2 subplots. Params
%                             might need to be set manually
% 
% -------------------------------------------------------------------------
%                           MISC
% ------------------------------------------------------------------------- 
% fixepsbbox                - fixes .eps bounding box to reduce extra whitespace
% suplabel                  - super label for subplots
% plotyyy                   - 3 y axes
% distinguishable_colors    - cycles best possible colors
% legendflex                - flexible legends
% scrollsubplot             - subplots with scrollbars.
% magnifyOnFigure           - is what is says
% pptfigure                 - 
% export_fig
%
% find_approx               - find approximate equality
% hist2d                    - plot 2d histogram
% matlab2latex              - is what it says
% run_avg2                  - Tom's running average code
% llpa                      - Ken's LLPA low pass filter for detiding.
% vgrid                     - Ken's ROMS vertical grid code.
% spice                     - calculate spiciness
% intersections             - find intersections of two curves
%
% fdep                      - lists dependencies for input .m file
% exportToZip
% v2struct
%
% 
% -------------------------------------------------------------------------
%                       TOOLBOXES
% -------------------------------------------------------------------------
% timeutil                  - time / date utilities 
%                             http://home.online.no/~pjacklam/matlab/software/util/index.html
% timeplt                   - 
% movieman                  - Ryan's movie stuff
% ROMS                      - roms_wilkin, arango, ROMS_tools
% tsplot                    - 
% cm_and_cb_utilities       - colormap and colorbar stuff

help dctools