%    [theta] = principalaxis(u,v)
%
% Returns principal axis angle in degrees calculated as
%    theta = 0.5*atand(2*cov(u,v)/(var(u) - var(v)))
%
% Generally you want to rotate this angle by -90 degrees so as to maximise
% variance or just use code below and figure which max you want
%
% angles = [theta 90+theta 180+theta 270+theta];
% for i = 1:length(angles)
%     vu(i) = var( up*cosd(angles(i)) - vp*sind(angles(i)));
%     vv(i) = var( up*sind(angles(i)) + vp*cosd(angles(i)));
% end

function [theta] = principalaxis(u,v)

    cuv   = cov(u,v);
    theta = 0.5*atand(2*cuv(2)/(var(u) - var(v)));