function [] = savepdf(filename)

    hfig = gcf;

    try
        hash = githash;
    catch ME
        hash = ' ';
        warning('githash didn''t work');
    end
    if strfind(hash, '(or any parent')
        hash = ' ';
    end

    % make later detection easier
    hash = ['hash:' hash];
    str = getAnnotation(gcf);
    hash = [hash ' | ' str];

    hfig.PaperPositionMode = 'auto';
    hfig.PaperSize = hfig.PaperPosition(3:4);
    print(hfig, '-dpdf', filename);
    system(['exiftool -overwrite_original -Producer=''' ...
            hash ''' ' filename]);