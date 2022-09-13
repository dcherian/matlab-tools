%       
%       function [freq,pow,cn] = dcfft(data,dt)
%              - Returns frequency, power, coefficients (NOT to be scaled by 2/N) returned by MATLAB's fft
%              - Takes time series and delta_t between 2 elements in the series
%              - folds over the engative frequencies
%              - cn = a_n + i*b_n

% modified from http://www.mathworks.com/support/tech-notes/1700/1702.html


function [freq,pow,cn] = dcfft(data,dt)

    % Sampling frequency 
    Fs = 1/dt; 
    x = data;    
    
    % Wunsch recommendation (time series primer page 38)
    if isprime(length(data))
        fprintf('\n Warning number of points is prime. Truncating by 1\n\n');
        x = x(1:end-1); 
    end

    % Use next highest power of 2 greater than or equal to length(x) to calculate FFT.
    nfft = length(x);%2^(nextpow2(length(x))); 

    % Take fft, padding with zeros so that length(fftx) is equal to nfft 
    fftx = fft(x);%,nfft); 

    % Calculate the numberof unique points
    NumUniquePts = ceil((nfft+1)/2); 

    % FFT is symmetric, throw away second half 
    fftx = fftx(1:NumUniquePts); 
    cn = 2/length(data) * fftx;

    % Take the magnitude of fft of x and scale the fft so that it is not a function of the length of x
    mx = abs(fftx)/length(x); 

    % Take the square of the magnitude of fft of x. 
    mx = mx.^2;   

    % Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
    % The DC component and Nyquist component (extra when N is even),
    % if it exists, are unique and should not be multiplied by 2.
    if rem(nfft, 2) % odd nfft excludes Nyquist point
      mx(2:end) = mx(2:end)*2;
    else
      mx(2:end-1) = mx(2:end-1)*2;
    end    
    pow = mx;
    
    % This is an evenly spaced frequency vector with NumUniquePts points. 
    freq = (0:NumUniquePts-1)*Fs/nfft;
    
    % Check whether parseval's theorem is satisfied
    %chkparsvl(data,cn,dt,1/length(data));
    
% Generate the plot, title and labels. 

% plot(f,mx); 
% title('Power Spectrum of a 200Hz Sine Wave'); 
% xlabel('Frequency (Hz)'); 
% ylabel('Power');
    
%      N = length(data);
%      
%      cn = fft(data);
%      cn = cn(1:N/2);
%      pow = (abs(cn)*2/N).^2;
%      
%      freq = [0:N/2-1]/(dt)/(N);
     
%     % Use next highest power of 2 greater than or equal to length(data) to calculate FFT.
%     % . For the fastest possible ffts, you will want to pad your data with enough zeros 
%     % to make its length a power of 2. The built-in FFT function does this for you automatically, 
%     % if you give a second argument specifying the overall length of the fft
%     nfft= 2^(nextpow2(length(data))); 
%    
%     cn = fft(data,nfft); 
%     
%     % FFT is symmetric, throw away second half
%     NumUniquePts = ceil((nfft+1)/2); 
%     cn = cn(1:NumUniquePts);
%     
%     pow = (abs(cn)*2/nfft).^2;
%     
% 
%     freq = (0:NumUniquePts-1)/dt;%/nfft
    
    
    % Don't understand this part
    % Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
    % The DC component and Nyquist component, if it exists, are unique and should not be multiplied by 2.
%     if rem(nfft, 2) % odd nfft excludes Nyquist point
%       pow(2:end) = pow(2:end)*2;
%     else
%       pow(2:end -1) = pow(2:end -1)*2;
%     end
    
    