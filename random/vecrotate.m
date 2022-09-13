% Rotates the vectors x and y into a new co-ordinate system
%       [rotx,roty] = vecrotate(x,y,theta)
%               theta - in degrees

function [rotx,roty] = vecrotate(x,y,theta)

    rotx =  x*cosd(theta) - y*sind(theta);
    roty =  x*sind(theta) + y*cosd(theta);