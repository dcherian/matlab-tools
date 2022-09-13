function LowerLines
    obj = findobj(gca, 'Tag', 'dcline');
    if ~isempty(obj)
        uistack(obj, 'bottom');
    end
end