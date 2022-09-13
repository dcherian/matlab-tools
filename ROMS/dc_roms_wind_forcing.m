function [] = dc_roms_wind_forcing(S,frcname)

    error('NEED TO MAKE THIS XI_V and XI_U etc');
    if exist(frcname,'file'), delete(frcname); end

    nccreate(frcname,'sms_time','Dimensions',{'sms_time' Inf});
    nccreate(frcname,'sustr','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2 'sms_time'});
    nccreate(frcname,'svstr','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2 'sms_time'});
    
    nccreate(frcname,'x_rho','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2'});
    nccreate(frcname,'y_rho','Dimensions',{'xi_rho' S.Lm+2 'eta_rho' S.Mm+2});
    
    ncwrite(frcname,'x_rho',S.x_rho);
    ncwrite(frcname,'y_rho',S.y_rho);
    
    % add attributes
    ncwriteatt(frcname,'sustr','long_name','surface u-momentum stress');
    ncwriteatt(frcname,'sustr','units', 'Newton meter-2');
    ncwriteatt(frcname,'sustr','coordinates', 'x_rho y_rho sms_time');
    ncwriteatt(frcname,'sustr','time','sms_time');
    
    ncwriteatt(frcname,'svstr','long_name','surface v-momentum stress');
    ncwriteatt(frcname,'svstr','units', 'Newton meter-2');
    ncwriteatt(frcname,'svstr','coordinates', 'x_rho y_rho sms_time');
    ncwriteatt(frcname,'svstr','time','sms_time');
    
    ncwriteatt(frcname,'sms_time','long_name','seconds since 0001-01-01 00:00:00');
    ncwriteatt(frcname,'sms_time','units','seconds');
    
    ncwriteatt(frcname,'x_rho','long_name','X-location of RHO-points');
    ncwriteatt(frcname,'x_rho','units','meter');
    ncwriteatt(frcname,'y_rho','long_name','Y-location of RHO-points');
    ncwriteatt(frcname,'y_rho','units','meter');
