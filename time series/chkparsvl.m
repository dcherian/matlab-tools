% Checks whether parseval's relation is satisfied
%       [res] = chkparsvl(data, data_fft, dt, ds)

function [] = chkparsvl(data, data_fft, dt, ds)

    data_sum = var(data);%trapz((data.^2))*dt/length(data);%
    fft_sum = sum(0.5*abs(data_fft).^2);%trapz(abs((data_fft)).^2)*ds;%
    
    if fft_sum ~= data_sum
        fprintf('\n\n Parsevals? Difference = %f = %f %% \n\n', (fft_sum-data_sum),(fft_sum-data_sum)/max(data_sum,fft_sum) * 100);
    end