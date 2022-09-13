% plot 'var1' at appropriate depths 'z_series' for x_series
%       [] = disp_plot(var1,x_series,z_series)

function [] = disp_plot(var1,x_series,z_series)
    
	if ndims(x_series)~=2 || ndims(z_series)~=2
        fprintf(' \n Error in either X/Z series. ');
        return;
	end
    
    if size(x_series,2) ~=1
        x_series = x_series';
    end
    if size(z_series,2) ~=1
        z_series = z_series';
    end
    
    figure
    
    % detrend i.e., remove mean
    m = nanmean(var1);
    var_dt = var1-repmat(m,length(var1),1);
	MM = repmat(z_series',size(var1,1),1);
	eta1 = var_dt + MM;
	
    plot(x_series,eta1);
	set(gca,'ydir','reverse');
    ylabel('Z(m)');