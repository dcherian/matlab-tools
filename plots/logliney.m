function [] = logliney(y,factor)
    
    hFig = evalin('caller','gcf');
    hAxis = evalin('caller','gca');
    
    figure(hFig);
    hold on;
    
    xax = (get(hAxis,'XLim'));
    
    %loglog([x; x]',repmat(yax,length(x),1),'k-','LineWidth',1.5);
    
    for i=1:length(y)       
        loglog(xax,[y(i) y(i)],'k--','LineWidth',1.5);
        %if mod(i,2)
            text(xax(1),y(i),num2str(y(i)),'Rotation',0,'VerticalAlignment','Bottom','FontWeight','Bold');
        %else
        %    text(x(i),yax(mod(i,2)+2),num2str(1./factor./x(i)),'Rotation',90,'VerticalAlignment','Bottom','FontWeight','Bold');
        %end
    end