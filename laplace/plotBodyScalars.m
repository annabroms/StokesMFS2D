function plotBodyScalars(q,R,field_left,field_right,font_size,plot_titles,cbar_labels,figure_title)
%PLOTBODYSCALARS Visualize two scalar fields on circular bodies.
%
%   plotBodyScalars(q,R,field_left,field_right)
%   plotBodyScalars(...,font_size)
%   plotBodyScalars(...,font_size,plot_titles,cbar_labels,figure_title)
%
% The default labels preserve the capacitance-style visualization used by
% solve_cap_peanut. Callers may override the labels to reuse the helper for
% other bodywise scalar fields.

q = q(:);
field_left = field_left(:);
field_right = field_right(:);
P = numel(q);

if nargin < 5 || isempty(font_size)
    font_size = 14;
end
if nargin < 6 || isempty(plot_titles)
    plot_titles = {'Given body voltages','Computed net charges'};
end
if nargin < 7 || isempty(cbar_labels)
    cbar_labels = {'Voltage','Net charge'};
end
if nargin < 8 || isempty(figure_title)
    figure_title = 'Bodywise voltages and net charges';
end

if isscalar(R)
    rad = repmat(R,P,1);
else
    rad = R(:);
end

theta = linspace(0,2*pi,200);
xmin = min(real(q)-rad);
xmax = max(real(q)+rad);
ymin = min(imag(q)-rad);
ymax = max(imag(q)+rad);
pad = 0.1*max(rad);

figure();
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

vals = {field_left, field_right};
titles = plot_titles;
cmaps = {parula(256),blueWhiteRedMap(256)};

for it = 1:2
    ax = nexttile;
    hold(ax,'on');

    for k = 1:P
        zk = q(k) + rad(k)*(cos(theta)+1i*sin(theta));
        fill(ax,real(zk),imag(zk),vals{it}(k), ...
            'EdgeColor',[0.2 0.2 0.2],'LineWidth',0.75);
    end

    axis(ax,'equal');
    xlim(ax,[xmin-pad xmax+pad]);
    ylim(ax,[ymin-pad ymax+pad]);
    xlabel(ax,'$x$','Interpreter','latex');
    ylabel(ax,'$y$','Interpreter','latex');
    title(ax,titles{it},'Interpreter','latex');
    box(ax,'on');
    colormap(ax,cmaps{it});
    ax.FontSize = font_size;
    ax.TitleFontSizeMultiplier = 1.0;
    ax.LabelFontSizeMultiplier = 1.0;
    ax.TickLabelInterpreter = 'latex';

    if min(vals{it}) < 0 && max(vals{it}) > 0
        vmax = max(abs(vals{it}));
        clim(ax,[-vmax vmax]);
    end

    c = colorbar(ax);
    c.Label.String = cbar_labels{it};
    c.FontSize = font_size;
    c.Label.FontSize = font_size;
    c.TickLabelInterpreter = 'latex';
    c.Label.Interpreter = 'latex';
end

sgtitle(figure_title,'FontSize',font_size,'Interpreter','latex');
end

function cmap = blueWhiteRedMap(n)
if nargin < 1 || isempty(n)
    n = 256;
end

x = linspace(0,1,n)';
cmap = zeros(n,3);

mid = 0.5;
left = x <= mid;
right = x > mid;

tleft = x(left)/mid;
tright = (x(right)-mid)/(1-mid);

cmap(left,1) = tleft;
cmap(left,2) = tleft;
cmap(left,3) = 0.3 + 0.7*tleft;

cmap(right,1) = 1.0;
cmap(right,2) = 1.0 - 0.8*tright;
cmap(right,3) = 1.0 - tright;
end
