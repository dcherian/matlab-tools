function [] = ColumnLegendFontSize(legobj, fs)

    objs = findobj(legobj, 'Type', 'Text');

    for ii=1:length(objs)
        objs(ii).FontSize = fs;
    end

end