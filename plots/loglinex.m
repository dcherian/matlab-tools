% Plots vertical line at a given x (can be a vector) and a factor (can also
% be vector) to scale the labelled 1./x value by.
% Eg: if 'x' is freq in cph and factor is 24, label shows PERIOD in days.
% does num2str(1./factor./x(i))
%       [] = loglinex(x, factor)

function [] = loglinex(x,factor)
    
    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    figure(hFig);
    hold on;
    
    yax = (get(hAxis,'YLim'));
    
    %loglog([x; x]',repmat(yax,length(x),1),'k-','LineWidth',1.5);
    
    for i=1:length(x)       
        loglog([x(i) x(i)],yax,'k--','LineWidth',1.5);
        %if mod(i,2)
            text(x(i),yax(1),num2str(1./factor./x(i)),'Rotation',90,'VerticalAlignment','Bottom','FontWeight','Bold');
        %else
        %    text(x(i),yax(mod(i,2)+2),num2str(1./factor./x(i)),'Rotation',90,'VerticalAlignment','Bottom','FontWeight','Bold');
        %end
    end