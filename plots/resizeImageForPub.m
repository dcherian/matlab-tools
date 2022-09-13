%  [] = resizeImageForPub(option)
%     option == 'test' || 'onecolumn' || 'potrait' || 'landscape'
function [] = resizeImageForPub(option)

    if strcmpi(option, 'test')
        test_resize;
    end

    hfig = gcf;
    funits = hfig.Units;

    hfig.Units = 'inches';
    hfig.Resize = 'off'; % for some reason required.

    if strcmpi(option, 'onecolumn')
        hfig.Position(3) = (8.5-1)/2;
    end

    if strcmpi(option, 'portrait')
        hfig.Position(3) = 8.5-1;
    end

    if strcmpi(option, 'onecolumnfullpage')
        hfig.Position(3) = (8.5-1)/2;
        hfig.Position(4) = 11.5-2;
        hfig.Units = 'normalized';
        hfig.Position(2) = 0.9 - hfig.Position(4);
    end

    if strcmpi(option, 'portraitfullpage')
        hfig.Position(3) = 8.5-1;
        hfig.Position(4) = 11.5-2;
        hfig.Units = 'normalized';
        hfig.Position(2) = 0.9 - hfig.Position(4);
    end

    if strcmpi(option, 'landscape')
        hfig.Position(3) = 11.5-1;
    end

    hfig.Units = funits;
    pause(1);
    hfig.Resize = 'on';
end

function test_resize

    figure; hfig = gcf;
    hfig.Units = 'inches';
    plot(rand([10 10]));
    pbaspect([1.617 1 1]);
    title('asdasd'); xlabel('X'); ylabel('Y');
    beautify;
    drawnow;
    resizeImageForPub('portrait');
    assert(abs(hfig.Position(3) - (8.5-1)) < 10*eps);
    close;
end