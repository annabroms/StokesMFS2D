clear;
close all;
clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
addpath(genpath(repoRoot));
set(0,'DefaultFigureVisible','off');

run(fullfile(repoRoot,'visuals','matlab_gen','mar26_visualse_coarse_fine.m'));

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(scriptDir));
annotationVarName = 'mar26VisualiseCoarseFineAnnotations';
ann = evalin('base',annotationVarName);

fig1Name = 'mar26 visualise coarse fine - coarse cluster';
fig2Name = 'mar26 visualise coarse fine - fine chain';
fig3Name = 'mar26 visualise coarse fine - closest neighbours';

pdfFile1 = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig1.pdf');
annotationTexFile1 = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig1_annotations.tex');
demoTexFile1 = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig1_demo.tex');
pdfFile2 = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig2.pdf');
annotationTexFile2 = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig2_annotations.tex');
demoTexFile2 = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig2_demo.tex');
pdfFile = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig3.pdf');
annotationTexFile = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig3_annotations.tex');
demoTexFile = fullfile(scriptDir,'mar26_visualse_coarse_fine_fig3_demo.tex');

fig1ExportScale = getAxisDataScale(fig1Name,[]);
exportCleanFullFigure(fig1Name,pdfFile1);
exportCleanFigure(fig2Name,pdfFile2,ann.secondFigureBounds,fig1ExportScale);
exportCleanFigure(fig3Name,pdfFile,[],[]);
writeEmptyAnnotationTex(annotationTexFile1);
writeDemoTex(demoTexFile1,'mar26_visualse_coarse_fine_fig1.pdf', ...
    'mar26_visualse_coarse_fine_fig1_annotations.tex');
writeFigure2AnnotationTex(annotationTexFile2,ann);
writeDemoTex(demoTexFile2,'mar26_visualse_coarse_fine_fig2.pdf', ...
    'mar26_visualse_coarse_fine_fig2_annotations.tex');
writeAnnotationTex(annotationTexFile,ann);
writeDemoTex(demoTexFile,'mar26_visualse_coarse_fine_fig3.pdf', ...
    'mar26_visualse_coarse_fine_fig3_annotations.tex');

fprintf('Wrote figure PDF: %s\n',pdfFile1);
fprintf('Wrote annotation overlay: %s\n',annotationTexFile1);
fprintf('Wrote standalone demo: %s\n',demoTexFile1);
fprintf('Wrote figure PDF: %s\n',pdfFile2);
fprintf('Wrote annotation overlay: %s\n',annotationTexFile2);
fprintf('Wrote standalone demo: %s\n',demoTexFile2);
fprintf('Wrote figure PDF: %s\n',pdfFile);
fprintf('Wrote annotation overlay: %s\n',annotationTexFile);
fprintf('Wrote standalone demo: %s\n',demoTexFile);


function exportCleanFigure(figName,pdfFile,bounds,targetScale)
[fig,ax] = findFigureAxes(figName);

if nargin >= 4 && ~isempty(bounds) && ~isempty(targetScale)
    setAxesPhysicalScale(fig,ax,bounds,targetScale);
end

delete(findall(ax,'Type','text'));
exportgraphics(ax,pdfFile,'ContentType','vector');
end


function exportCleanFullFigure(figName,pdfFile)
fig = findFigureByName(figName);
delete(findall(fig,'Type','text'));
exportgraphics(fig,pdfFile,'ContentType','vector');
end


function scale = getAxisDataScale(figName,bounds)
[fig,ax] = findFigureAxes(figName);
if nargin < 2 || isempty(bounds)
    bounds = [xlim(ax), ylim(ax)];
end
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


function [figPos,axPos] = getFigureAndAxesPositionInches(fig,ax)
figUnits = get(fig,'Units');
axUnits = get(ax,'Units');
cleanupObj = onCleanup(@() restoreUnits(fig,ax,figUnits,axUnits));

set(fig,'Units','inches');
set(ax,'Units','inches');
figPos = get(fig,'Position');
axPos = get(ax,'Position');
end


function restoreUnits(fig,ax,figUnits,axUnits)
set(fig,'Units',figUnits);
set(ax,'Units',axUnits);
end


function [fig,ax] = findFigureAxes(figName)
fig = findFigureByName(figName);
ax = findall(fig,'Type','axes');
if isempty(ax)
    error('Could not find axes in figure "%s".',figName);
end
ax = ax(1);
end


function fig = findFigureByName(figName)
fig = findall(groot,'Type','figure','Name',figName);
if isempty(fig)
    error('Could not find figure "%s".',figName);
end
fig = fig(1);
end


function writeFigure2AnnotationTex(filename,ann)
bounds = ann.secondFigureBounds;
xmin = bounds(1);
xmax = bounds(2);
ymin = bounds(3);
ymax = bounds(4);

fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Auto-generated from mar26_visualse_coarse_fine.m\n');
fprintf(fid,'\\begin{scope}[shift={(img.south west)},\n');
fprintf(fid,'  x={($(img.south east)-(img.south west)$)},\n');
fprintf(fid,'  y={($(img.north west)-(img.south west)$)}]\n');
writeNodeGroup(fid,ann.figure2CircleNumberLabels,'mar26 number',xmin,xmax,ymin,ymax,false,false);
fprintf(fid,'\\end{scope}\n');
end


function writeEmptyAnnotationTex(filename)
fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end

cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%% Auto-generated from mar26_visualse_coarse_fine.m\n');
end


function writeAnnotationTex(filename,ann)
bounds = ann.thirdFigureBounds;
xmin = bounds(1);
xmax = bounds(2);
ymin = bounds(3);
ymax = bounds(4);

fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end

cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid,'%% Auto-generated from mar26_visualse_coarse_fine.m\n');
fprintf(fid,'%% Data bounds: xmin=%.16g, xmax=%.16g, ymin=%.16g, ymax=%.16g\n',xmin,xmax,ymin,ymax);
fprintf(fid,'%% Normalization used here:\n');
fprintf(fid,'%%   u = (x - xmin)/(xmax - xmin)\n');
fprintf(fid,'%%   v = (y - ymin)/(ymax - ymin)\n');
fprintf(fid,'\\begin{scope}[shift={(img.south west)},\n');
fprintf(fid,'  x={($(img.south east)-(img.south west)$)},\n');
fprintf(fid,'  y={($(img.north west)-(img.south west)$)}]\n');

writeNodeGroup(fid,ann.psiLabels,'mar26 psi',xmin,xmax,ymin,ymax,true,false);
writeNodeGroup(fid,ann.phiLabels,'mar26 phi',xmin,xmax,ymin,ymax,true,false);
writeNodeGroup(fid,ann.topPlusLabels,'mar26 top plus',xmin,xmax,ymin,ymax,true,false);
writeNodeGroup(fid,ann.verticalPlusLabels,'mar26 vertical plus',xmin,xmax,ymin,ymax,true,false);
writeNodeGroup(fid,ann.equalsLabels,'mar26 equals',xmin,xmax,ymin,ymax,true,true);
writeNodeGroup(fid,ann.pairTitleLabels,'mar26 pair title',xmin,xmax,ymin,ymax,false,false);
writeNodeGroup(fid,ann.pairFocusLabels,'mar26 pair focus',xmin,xmax,ymin,ymax,false,false);
writeNodeGroup(fid,ann.pairNeighbourLabels,'mar26 pair neighbour',xmin,xmax,ymin,ymax,false,false);
writeNodeGroup(fid,ann.circleNumberLabels,'mar26 number',xmin,xmax,ymin,ymax,false,false);
writeNodeGroup(fid,ann.neighbourCircleNumberLabels,'mar26 neighbour number',xmin,xmax,ymin,ymax,false,false);

fprintf(fid,'\\end{scope}\n');
end


function writeNodeGroup(fid,nodes,styleName,xmin,xmax,ymin,ymax,useMathMode,rotateEquals)
for k = 1:numel(nodes)
    pos = nodes(k).position;
    u = (pos(1) - xmin)/(xmax - xmin);
    v = (pos(2) - ymin)/(ymax - ymin);
    nodeText = nodes(k).text;

    if useMathMode
        content = sprintf('$%s$',nodeText);
    else
        content = nodeText;
    end

    if rotateEquals
        fprintf(fid,'  \\node[%s,rotate=90] at (%.16f, %.16f) {%s};\n', ...
            styleName,u,v,content);
    else
        fprintf(fid,'  \\node[%s] at (%.16f, %.16f) {%s};\n', ...
            styleName,u,v,content);
    end
end
end


function writeDemoTex(filename,pdfName,annotationName)
fid = fopen(filename,'w');
if fid < 0
    error('Could not open %s for writing.',filename);
end

cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid,'\\documentclass[tikz,border=4pt]{standalone}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage{tikz}\n');
fprintf(fid,'\\usetikzlibrary{calc}\n');
fprintf(fid,'\n');
fprintf(fid,'%% Adjust these font commands later without touching the saved positions.\n');
fprintf(fid,'\\newcommand{\\marTwentySixFormulaFont}{\\fontsize{14}{14}\\selectfont}\n');
fprintf(fid,'\\newcommand{\\marTwentySixNumberFont}{\\fontsize{12}{12}\\selectfont}\n');
fprintf(fid,'\n');
fprintf(fid,'\\tikzset{\n');
fprintf(fid,'  mar26 psi/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 phi/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 top plus/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 vertical plus/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 equals/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 pair title/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 pair focus/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 pair neighbour/.style={font=\\marTwentySixFormulaFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt},\n');
fprintf(fid,'  mar26 number/.style={font=\\marTwentySixNumberFont,text=black,anchor=center,inner sep=0pt,outer sep=0pt,text height=1.5ex,text depth=0.25ex},\n');
fprintf(fid,'  mar26 neighbour number/.style={font=\\marTwentySixNumberFont,text={rgb,1:red,0.55;green,0.55;blue,0.55},anchor=center,inner sep=0pt,outer sep=0pt,text height=1.5ex,text depth=0.25ex}\n');
fprintf(fid,'}\n');
fprintf(fid,'\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\begin{tikzpicture}\n');
fprintf(fid,'  \\node[inner sep=0pt,outer sep=0pt] (img) at (0,0) {%%\n');
fprintf(fid,'    \\includegraphics{%s}%%\n',pdfName);
fprintf(fid,'  };\n');
fprintf(fid,'  \\input{%s}\n',annotationName);
fprintf(fid,'\\end{tikzpicture}\n');
fprintf(fid,'\\end{document}\n');
end
