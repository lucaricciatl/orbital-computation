function [fig,ax] = set_graphichs_properties()
    fig = figure;
    ax=axes;
	opengl software
    hold on
    light
    axis tight
    shading interp
    lighting gouraud 
    grid on
    xlabel('x[au]')
    ylabel('y[au]')
    zlabel('z[au]')
    set(gca,'Color','k')
    set(gcf,'Color','black');
    set(gca,'xcolor','w')
    set(gca,'ycolor','w')
    set(gca,'zcolor','w')
    view(3)
    width=1920/2;
    height=1080/2;
    set(gcf,'position',[0,0,width,height])
    xlim([-11 11])
    ylim([-11 11])
    zlim([-2 2])
    daspect([1 1 1])
    set(groot,'DefaultFigureGraphicsSmoothing','on')

end