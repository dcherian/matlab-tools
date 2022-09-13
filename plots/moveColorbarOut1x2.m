% provide handle to last colorbar in 1x2 suplot grid
function [] = moveColorbar1x2(hcb)
    axpos = get(gca, 'Position');
    hcb.Position(1) = 0.92;
    hcb.Position(2) = 0.35;
    hcb.Position(3) = 0.017;
    hcb.Position(4) = axpos(4)/3;
end