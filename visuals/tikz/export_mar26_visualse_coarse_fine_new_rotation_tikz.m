clear;
close all;
clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(genpath(repoRoot));

referenceScales = getUnrotatedReferenceScales(repoRoot);
out = mar26_visualse_coarse_fine_new_rotation();

pdfFile4 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig4.pdf');
annotationTexFile4 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig4_annotations.tex');
demoTexFile4 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig4_demo.tex');

pdfFile5 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig5.pdf');
annotationTexFile5 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig5_annotations.tex');
demoTexFile5 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig5_demo.tex');

pdfFile6 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig6.pdf');
annotationTexFile6 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig6_annotations.tex');
demoTexFile6 = fullfile(scriptDir,'mar26_visualse_coarse_fine_new_rotation_fig6_demo.tex');

fig4Nodes = collectStyledText(out.figure4Axes);
fig4Nodes = promoteMathBold(fig4Nodes);
fig4Bounds = [xlim(out.figure4Axes), ylim(out.figure4Axes)];
fig6Nodes = collectStyledText(out.figure6Axes);
fig6Nodes = promoteMathBold(fig6Nodes);
fig6Bounds = [xlim(out.figure6Axes), ylim(out.figure6Axes)];

exportCleanAxesFigure(out.figure4,out.figure4Axes,pdfFile4,fig4Bounds, ...
    referenceScales.closestNeighboursScale);
exportScaledFigureWithInset(out.figure5,out.figure5MainAxes,out.figure5MainBounds, ...
    out.figure5InsetAxes,out.figure5InsetBounds,out.figure5InsetGapInches, ...
    out.figure5InsetCenterOffsetInches,out.figure5InsetBottomMarginInches, ...
    pdfFile5,referenceScales.coarseClusterScale);
exportCleanAxesFigure(out.figure6,out.figure6Axes,pdfFile6,fig6Bounds, ...
    referenceScales.closestNeighboursScale);

writeAnnotationTex(annotationTexFile4,fig4Nodes,fig4Bounds);
writeEmptyAnnotationTex(annotationTexFile5);
writeAnnotationTex(annotationTexFile6,fig6Nodes,fig6Bounds);
writeDemoTex(demoTexFile4,'mar26_visualse_coarse_fine_new_rotation_fig4.pdf', ...
    'mar26_visualse_coarse_fine_new_rotation_fig4_annotations.tex');
writeDemoTex(demoTexFile5,'mar26_visualse_coarse_fine_new_rotation_fig5.pdf', ...
    'mar26_visualse_coarse_fine_new_rotation_fig5_annotations.tex');
writeDemoTex(demoTexFile6,'mar26_visualse_coarse_fine_new_rotation_fig6.pdf', ...
    'mar26_visualse_coarse_fine_new_rotation_fig6_annotations.tex');

fprintf('Wrote figure PDF: %s\n',pdfFile4);
fprintf('Wrote annotation overlay: %s\n',annotationTexFile4);
fprintf('Wrote standalone demo: %s\n',demoTexFile4);
fprintf('Wrote figure PDF: %s\n',pdfFile5);
fprintf('Wrote annotation overlay: %s\n',annotationTexFile5);
fprintf('Wrote standalone demo: %s\n',demoTexFile5);
fprintf('Wrote figure PDF: %s\n',pdfFile6);
fprintf('Wrote annotation overlay: %s\n',annotationTexFile6);
fprintf('Wrote standalone demo: %s\n',demoTexFile6);


function nodes = collectStyledText(ax)
textHandles = findall(ax,'Type','text');
nodes = struct('text',{},'position',{},'style',{});

for k = 1:numel(textHandles)
    textNow = textHandles(k);
    textString = get(textNow,'String');
    textPos = get(textNow,'Position');
    textColor = get(textNow,'Color');

    node = struct();
    node.text = textString;
    node.position = textPos(1:2);
    node.style = classifyTextStyle(textString,textColor);
    nodes(end+1) = node; %#ok<AGROW>
end

nodes = flip(nodes);
end


function nodes = promoteMathBold(nodes)
for k = 1:numel(nodes)
    switch nodes(k).style
        case 'mar26 psi'
            nodes(k).text = strrep(nodes(k).text,'\psi','\boldsymbol{\psi}');
        case 'mar26 phi'
            nodes(k).text = strrep(nodes(k).text,'\phi','\boldsymbol{\phi}');
        case 'mar26 pair title'
            nodes(k).text = strrep(nodes(k).text,'\chi','\boldsymbol{\chi}');
        case 'mar26 pair data'
            nodes(k).text = strrep(nodes(k).text,'-\phi','-\boldsymbol{\phi}');
    end
end
end


function exportCleanAxesFigure(fig,ax,pdfFile,bounds,targetScale)
delete(findall(fig,'Type','text'));
drawnow;

if nargin >= 5 && ~isempty(bounds) && ~isempty(targetScale)
    setAxesPhysicalScale(fig,ax,bounds,targetScale);
end

exportgraphics(ax,pdfFile,'ContentType','vector');
end


function exportScaledFigure(fig,ax,pdfFile,bounds,targetScale)
if nargin >= 5 && ~isempty(bounds) && ~isempty(targetScale)
    setAxesPhysicalScale(fig,ax,bounds,targetScale);
end

drawnow;
exportgraphics(fig,pdfFile,'ContentType','vector');
end


function exportScaledFigureWithInset(fig,mainAx,mainBounds,insetAx,insetBounds, ...
    insetGapInches,insetCenterOffsetInches,insetBottomMarginInches,pdfFile,targetScale)
if nargin >= 10 && ~isempty(mainBounds) && ~isempty(targetScale)
    setAxesPhysicalScale(fig,mainAx,mainBounds,targetScale);
end

arrangeInsetAxesToMatchMainScale(fig,mainAx,mainBounds,insetAx,insetBounds, ...
    insetGapInches,insetCenterOffsetInches,insetBottomMarginInches);

drawnow;
exportgraphics(fig,pdfFile,'ContentType','vector');
end


function scales = getUnrotatedReferenceScales(repoRoot)
run(fullfile(repoRoot,'visuals','matlab_gen','mar26_visualse_coarse_fine.m'));
[clusterFig,clusterAx] = findFigureAxes('mar26 visualise coarse fine - coarse cluster');
clusterBounds = [xlim(clusterAx), ylim(clusterAx)];
[refFig,refAx] = findFigureAxes('mar26 visualise coarse fine - closest neighbours');
refBounds = [xlim(refAx), ylim(refAx)];
scales = struct();
scales.coarseClusterScale = getAxisDataScale(clusterFig,clusterAx,clusterBounds);
scales.closestNeighboursScale = getAxisDataScale(refFig,refAx,refBounds);
closeReferenceFigures();
end


function closeReferenceFigures()
figs = findall(groot,'Type','figure');
for k = 1:numel(figs)
    figName = get(figs(k),'Name');
    if startsWith(figName,'mar26 visualise coarse fine -') && ~contains(figName,'new rotation')
        close(figs(k));
    end
end
end


function scale = getAxisDataScale(fig,ax,bounds)
[~,axPos] = getFigureAndAxesPositionInches(fig,ax);
dx = bounds(2) - bounds(1);
dy = bounds(4) - bounds(3);
scale = min(axPos(3)/dx, axPos(4)/dy);
end


function setAxesPhysicalScale(fig,ax,bounds,targetScale)
[figPos,axPos] = getFigureAndAxesPositionInches(fig,ax);

marginLeft = axPos(1);
marginBottom = axPos(2);
marginRight = figPos(3) - axPos(1) - axPos(3);
marginTop = figPos(4) - axPos(2) - axPos(4);

desiredAxWidth = targetScale*(bounds(2) - bounds(1));
desiredAxHeight = targetScale*(bounds(4) - bounds(3));

figPos(3) = marginLeft + desiredAxWidth + marginRight;
figPos(4) = marginBottom + desiredAxHeight + marginTop;
axPos(3) = desiredAxWidth;
axPos(4) = desiredAxHeight;

set(fig,'Units','inches','Position',figPos);
set(ax,'Units','inches','Position',axPos);
drawnow;
end


function arrangeInsetAxesToMatchMainScale(fig,axMain,mainBounds,axInset,insetBounds, ...
    gapInches,centerOffsetInches,bottomMarginInches)
[figPos,mainPos,insetPos,figUnits,mainUnits,insetUnits] = ...
    getFigureAndAxesPositionsInches(fig,axMain,axInset);
cleanupObj = onCleanup(@() restoreFigureAxesUnits(fig,axMain,axInset, ...
    figUnits,mainUnits,insetUnits)); %#ok<NASGU>

mainScale = min(mainPos(3)/(mainBounds(2) - mainBounds(1)), ...
    mainPos(4)/(mainBounds(4) - mainBounds(3)));
insetWidth = mainScale*(insetBounds(2) - insetBounds(1));
insetHeight = mainScale*(insetBounds(4) - insetBounds(3));

insetX = mainPos(1) + 0.5*(mainPos(3) - insetWidth) + centerOffsetInches;
insetY = mainPos(2) - gapInches - insetHeight;

if insetY < bottomMarginInches
    shiftUp = bottomMarginInches - insetY;
    figPos(4) = figPos(4) + shiftUp;
    mainPos(2) = mainPos(2) + shiftUp;
    insetPos(2) = insetPos(2) + shiftUp;
    set(fig,'Units','inches','Position',figPos);
    set(axMain,'Units','inches','Position',mainPos);
    insetY = bottomMarginInches;
end

set(axInset,'Units','inches','Position',[insetX, insetY, insetWidth, insetHeight]);
drawnow;
end


function [figPos,axPos] = getFigureAndAxesPositionInches(fig,ax)
figUnits = get(fig,'Units');
axUnits = get(ax,'Units');
cleanupObj = onCleanup(@() restoreUnits(fig,ax,figUnits,axUnits)); %#ok<NASGU>

set(fig,'Units','inches');
set(ax,'Units','inches');
figPos = get(fig,'Position');
axPos = get(ax,'Position');
end


function [figPos,mainPos,insetPos,figUnits,mainUnits,insetUnits] = ...
    getFigureAndAxesPositionsInches(fig,axMain,axInset)
figUnits = get(fig,'Units');
mainUnits = get(axMain,'Units');
insetUnits = get(axInset,'Units');

set(fig,'Units','inches');
set(axMain,'Units','inches');
set(axInset,'Units','inches');

figPos = get(fig,'Position');
mainPos = get(axMain,'Position');
insetPos = get(axInset,'Position');
end


function restoreUnits(fig,ax,figUnits,axUnits)
set(fig,'Units',figUnits);
set(ax,'Units',axUnits);
end


function restoreFigureAxesUnits(fig,axMain,axInset,figUnits,mainUnits,insetUnits)
set(fig,'Units',figUnits);
set(axMain,'Units',mainUnits);
set(axInset,'Units',insetUnits);
end


function [fig,ax] = findFigureAxes(figName)
fig = findall(groot,'Type','figure','Name',figName);
if isempty(fig)
    error('Could not find figure "%s".',figName);
end
fig = fig(1);

ax = findall(fig,'Type','axes');
if isempty(ax)
    error('Could not find axes in figure "%s".',figName);
end
ax = ax(1);
end


function styleName = classifyTextStyle(textString,textColor)
if contains(textString,'\psi')
    styleName = 'mar26 psi';
elseif contains(textString,'\mathbf{u}')
    styleName = 'mar26 pair data';
elseif contains(textString,'\phi')
    styleName = 'mar26 phi';
elseif contains(textString,'\chi')
    styleName = 'mar26 pair title';
elseif strcmp(textString,'$+$')
    styleName = 'mar26 plus';
elseif strcmp(textString,'$=$')
    styleName = 'mar26 equals';
elseif all(isstrprop(textString,'digit'))
    if norm(textColor - [0.55 0.55 0.55]) < 1e-6
        styleName = 'mar26 neighbour number';
    else
        styleName = 'mar26 number';
    end
else
    styleName = 'mar26 formula';
end
end


function writeAnnotationTex(filename,nodes,bounds)
xmin = bounds(1);
xmax = bounds(2);
ymin = bounds(3);
ymax = bounds(4);

fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'%% Auto-generated from mar26_visualse_coarse_fine_new_rotation.m\n');
fprintf(fid,'\\begin{scope}[shift={(img.south west)},\n');
fprintf(fid,'  x={($(img.south east)-(img.south west)$)},\n');
fprintf(fid,'  y={($(img.north west)-(img.south west)$)}]\n');

for k = 1:numel(nodes)
    u = (nodes(k).position(1) - xmin)/(xmax - xmin);
    v = (nodes(k).position(2) - ymin)/(ymax - ymin);
    fprintf(fid,'  \\node[%s] at (%.16f, %.16f) {%s};\n', ...
        nodes(k).style,u,v,nodes(k).text);
end

fprintf(fid,'\\end{scope}\n');
end


function writeEmptyAnnotationTex(filename)
fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Auto-generated from mar26_visualse_coarse_fine_new_rotation.m\n');
end


function writeDemoTex(filename,pdfName,annotationName)
fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid,'\\documentclass[tikz,border=4pt]{standalone}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage{tikz}\n');
fprintf(fid,'\\usetikzlibrary{calc}\n\n');
fprintf(fid,'\\newcommand{\\marTwentySixFormulaFont}{\\fontsize{14}{14}\\selectfont}\n');
fprintf(fid,'\\newcommand{\\marTwentySixNumberFont}{\\fontsize{12}{12}\\selectfont}\n\n');
fprintf(fid,'\\tikzset{\n');
fprintf(fid,'  mar26 formula/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 psi/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 phi/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 plus/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 equals/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 pair title/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 pair data/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 number/.style={font=\\marTwentySixNumberFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt,text height=1.5ex,text depth=0.25ex},\n');
fprintf(fid,'  mar26 neighbour number/.style={font=\\marTwentySixNumberFont,text={rgb,1:red,0.55;green,0.55;blue,0.55},anchor=center,inner sep=0pt,outer sep=0pt,text height=1.5ex,text depth=0.25ex}\n');
fprintf(fid,'}\n\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{tikzpicture}\n');
fprintf(fid,'  \\node[inner sep=0pt,outer sep=0pt] (img) at (0,0) {%%\n');
fprintf(fid,'    \\includegraphics{%s}%%\n',pdfName);
fprintf(fid,'  };\n');
fprintf(fid,'  \\input{%s}\n',annotationName);
fprintf(fid,'\\end{tikzpicture}\n');
fprintf(fid,'\\end{document}\n');
end
