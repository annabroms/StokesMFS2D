function styleLinePlot(ax,font_size,xlab,ylab,ttl)
%STYLELINEPLOT Apply consistent LaTeX styling to a simple line plot.

if nargin < 1 || isempty(ax)
    ax = gca;
end

set(ax,'FontSize',font_size,'TickLabelInterpreter','latex');
xlabel(ax,xlab,'Interpreter','latex','FontSize',font_size);
ylabel(ax,ylab,'Interpreter','latex','FontSize',font_size);
title(ax,ttl,'Interpreter','latex','FontSize',font_size);
lgd = legend(ax,'Location','best');
set(lgd,'Interpreter','latex','FontSize',font_size);
grid(ax,'on');
axis(ax,'tight');
end
