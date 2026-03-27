clear;
close all;
clc;

% TikZ export workflow:
% 1. Run this script to build figures 1-3 and populate the
%    mar26VisualiseCoarseFineAnnotations struct in the base workspace.
% 2. Run visuals/tikz/export_mar26_visualse_coarse_fine_fig3_tikz.m to
%    export a clean PDF for figure 3, regenerate the TikZ overlay, and now
%    also export figure 2 for side-by-side LaTeX use.
% 3. Compile the generated demo files in visuals/tikz/ with pdflatex, for
%    example:
%       pdflatex -interaction=nonstopmode -halt-on-error ...
%           mar26_visualse_coarse_fine_fig2_demo.tex
%       pdflatex -interaction=nonstopmode -halt-on-error ...
%           mar26_visualse_coarse_fine_fig3_demo.tex
% 4. Toggle rotate_group below when you want the figure-2/figure-3 geometry
%    rotated by 90 degrees before exporting to TikZ.

repo_root = fileparts(fileparts(mfilename('fullpath')));
if ~isempty(repo_root)
    addpath(genpath(repo_root));
end

fprintf('=== Coarse/Fine Visualisation (Mar 26, 2026) ===\n');

rad = 1;
domainPadding = 1.8;
nCirclePts = 240;

% Change this once to recolor all fluid regions in the script.
fluidColor = [0.82 0.91 0.98];
bodyColor = [1 1 1];
bodyEdgeColor = [0.72 0.77 0.82];

coarseNodeCount = 10;
coarseNodeColor = [0.10 0.34 0.82];
coarseNodeSize = 44;

fineNodeCount = 30;
fineNodeColor = [0.83 0.21 0.18];
fineNodeSize = 12;
fadedFineNodeColor = blendWithWhite(fineNodeColor,0.68);
formulaFontSize = 14;
numberFontSize = 12;
annotationTextColor = [0 0 0];
neighbourNumberColor = [0.55 0.55 0.55];
numberFontName = 'Courier';
numberBaseOffset = [0.00, 0.00];
numberSixOffset = [-0.01, 0.00];
numberSevenOffset = [0.00, 0.00];
numberPosition = @(center,label) getCircleNumberPosition( ...
    center,label,numberBaseOffset,numberSixOffset,numberSevenOffset);
drawNumber = @(ax,center,label,color) drawCircleNumber( ...
    ax,center,label,color,numberFontSize,numberFontName, ...
    numberBaseOffset,numberSixOffset,numberSevenOffset);
pairAnnotationColor = fineNodeColor;
row2LabelYOffset = 1.85;
pairLabelYOffset = 2.15;
pairLabelXOffset = 0.70;
pairTitleYOffset = 0.80;
bottomPairTitleLift = 0.30;
bottomVerticalPlusLift = 0.35;
rotate_group = true;
showPairVelocityLabels = false;
drawAnnotations = true;
exportAnnotationData = true;
annotationVarName = 'mar26VisualiseCoarseFineAnnotations';

annotationData = struct();
annotationData.annotationVarName = annotationVarName;
annotationData.psiLabels = struct('text',{},'position',{});
annotationData.phiLabels = struct('text',{},'position',{});
annotationData.topPlusLabels = struct('text',{},'position',{});
annotationData.verticalPlusLabels = struct('text',{},'position',{});
annotationData.equalsLabels = struct('text',{},'position',{});
annotationData.pairFocusLabels = struct('text',{},'position',{});
annotationData.pairNeighbourLabels = struct('text',{},'position',{});
annotationData.pairTitleLabels = struct('text',{},'position',{});
annotationData.circleNumberLabels = struct('text',{},'position',{});
annotationData.neighbourCircleNumberLabels = struct('text',{},'position',{});
annotationData.figure2CircleNumberLabels = struct('text',{},'position',{});

%% Coarse cluster: exterior fluid and single-body interior
delta = 1;
d = rad*(2 + delta);
a1 = d;
a2 = d*(0.5 + 1i*sqrt(3)/2);

%q_cluster = [0; a1; -a1; a2; -a2; a1-a2; -a1+a2; a1+a2];
%q_cluster = q_cluster - mean(q_cluster);
P = 8;
q_cluster = grow_cluster(P,delta,2,1,[],false,false);
focusIdx = 1;

figure('Color','w','Name','mar26 visualise coarse fine - coarse cluster');

ax1 = subplot(1,2,1);
setupFluidAxes(ax1,q_cluster,rad,fluidColor,domainPadding);
hold(ax1,'on');
for k = 1:numel(q_cluster)
    drawDisk(ax1,q_cluster(k),rad,bodyColor,'none',1.0,nCirclePts);
    drawBoundaryNodes(ax1,q_cluster(k),rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
end


ax2 = subplot(1,2,2);
setupFluidAxes(ax2,q_cluster,rad,fluidColor,domainPadding);
hold(ax2,'on');
for k = 1:numel(q_cluster)
    if k == focusIdx
        fillColor = bodyColor;
    else
        fillColor = fluidColor;
    end
    drawDisk(ax2,q_cluster(k),rad,fillColor,'none',1.0,nCirclePts);
    drawBoundaryNodes(ax2,q_cluster(k),rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
end


%% Fine chain: five close bodies and two farther bodies
closeGap = 0.15;
dClose = rad*(2 + closeGap);
turnAngles = (pi/180)*[15, 70, -20, 55];

q_close = zeros(5,1);
for k = 2:numel(q_close)
    q_close(k) = q_close(k-1) + dClose*exp(1i*turnAngles(k-1));
end

q_far = [q_close(2) - (3.7 - 2.9i); q_close(5) + (4.0 - 1.1i)];
shiftCenter = mean(q_close);
q_close = q_close - shiftCenter;
q_far = q_far - shiftCenter;
if rotate_group
    q_close = 1i*q_close;
    q_far = 1i*q_far;
end
q_fine = [q_close; q_far];

figure('Color','w','Name','mar26 visualise coarse fine - fine chain');
ax3 = axes;
setupFluidAxes(ax3,q_fine,rad,fluidColor,domainPadding);
annotationData.secondFigureBounds = getBounds(q_fine,rad,domainPadding);
hold(ax3,'on');
for k = 1:numel(q_fine)
    drawDisk(ax3,q_fine(k),rad,bodyColor,'none',1.0,nCirclePts);
end
for k = 1:numel(q_close)
    drawBoundaryNodes(ax3,q_close(k),rad,fineNodeCount,fineNodeColor,fineNodeSize);
end
for k = 1:numel(q_far)
    drawBoundaryNodes(ax3,q_far(k),rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
end
figure2LabelPositions = [q_close(:); q_far(:)];
figure2LabelTexts = [{'2','3','4','5','6'},{'1','7'}];
for k = 1:numel(figure2LabelPositions)
    annotationData.figure2CircleNumberLabels(k).text = figure2LabelTexts{k};
    annotationData.figure2CircleNumberLabels(k).position = ...
        numberPosition(figure2LabelPositions(k),figure2LabelTexts{k});
    if drawAnnotations
        drawNumber(ax3,figure2LabelPositions(k),figure2LabelTexts{k},annotationTextColor);
    end
end

%% Closest-neighbour views along the close chain
panelSpacing = 6.5;
topRowY = 0;
psiLabelY = 4.2;
equalsY = -4.1;
coarseFocusRowY = -7.4;
phiLabelY = coarseFocusRowY + row2LabelYOffset;
%pairTopRowY = -14.9;
pairTopRowY = -13.9;
verticalPlusTopY = 0.45*(coarseFocusRowY + pairTopRowY);
pairBottomRowY = -23.2;
pairBottomRowY = -21.4;
verticalPlusBottomY = 0.49*(pairTopRowY + pairBottomRowY) + bottomVerticalPlusLift;
q_strip = zeros(size(q_close));
for k = 1:numel(q_close)
    q_strip(k) = (k-1)*panelSpacing;
end

topFarDisplays = [q_strip(1) - panelSpacing + 1i*topRowY; ...
    q_strip(end) + panelSpacing + 1i*topRowY];
topDisplayCenters = [topFarDisplays(1); q_strip + 1i*topRowY; topFarDisplays(2)];
row2DisplayCenters = topDisplayCenters + 1i*coarseFocusRowY;

allPanelCenters = [topDisplayCenters(:); row2DisplayCenters(:); ...
    topDisplayCenters + 1i*(psiLabelY - topRowY); ...
    topDisplayCenters + 1i*(equalsY - topRowY); ...
    row2DisplayCenters + 1i*row2LabelYOffset];
neighbourLists = cell(numel(q_close),1);
pairDisplayLists = cell(numel(q_close),1);
for k = 1:numel(q_close)
    neighbourLists{k} = closestNeighbours(q_close,k);
    pairDisplayLists{k} = neighbourLists{k};
    if ismember(k,[2 3 4]) && numel(pairDisplayLists{k}) == 2
        pairDisplayLists{k} = pairDisplayLists{k}([2 1]);
    end
end

topPairTitleY = -inf;
bottomPairTitleY = -inf;
for k = 1:numel(q_close)
    if ~isempty(pairDisplayLists{k})
        relTop = q_close(pairDisplayLists{k}(1)) - q_close(k);
        pairTopMid = q_strip(k) + 1i*pairTopRowY;
        focusTop = pairTopMid - 0.5*relTop;
        neighbourTop = pairTopMid + 0.5*relTop;
        topPairTitleY = max(topPairTitleY, ...
            max(imag([focusTop, neighbourTop])) + rad + pairTitleYOffset);
    end

    if numel(pairDisplayLists{k}) >= 2
        relBottom = q_close(pairDisplayLists{k}(2)) - q_close(k);
        pairBottomMid = q_strip(k) + 1i*pairBottomRowY;
        focusBottom = pairBottomMid - 0.5*relBottom;
        neighbourBottom = pairBottomMid + 0.5*relBottom;
        bottomPairTitleY = max(bottomPairTitleY, ...
            max(imag([focusBottom, neighbourBottom])) + rad + pairTitleYOffset + bottomPairTitleLift);
    end
end

for k = 1:numel(q_close)

    centerNow = q_strip(k) + 1i*topRowY;
    localCenters = q_close(neighbourLists{k}) - q_close(k) + centerNow;
    allPanelCenters = [allPanelCenters; centerNow; localCenters(:); ...
        q_strip(k) + 1i*verticalPlusTopY]; %#ok<AGROW>

    if ~isempty(pairDisplayLists{k})
        relTop = q_close(pairDisplayLists{k}(1)) - q_close(k);
        pairTopMid = q_strip(k) + 1i*pairTopRowY;
        focusTop = pairTopMid - 0.5*relTop;
        neighbourTop = pairTopMid + 0.5*relTop;
        allPanelCenters = [allPanelCenters; focusTop; neighbourTop; ...
            q_strip(k) + 1i*pairTopRowY; ...
            real(pairTopMid) + 1i*topPairTitleY]; %#ok<AGROW>
        if showPairVelocityLabels
            allPanelCenters = [allPanelCenters; ...
                focusTop - pairLabelXOffset - 1i*pairLabelYOffset; ...
                neighbourTop + pairLabelXOffset - 1i*pairLabelYOffset]; %#ok<AGROW>
        end
    end

    if numel(pairDisplayLists{k}) >= 2
        relBottom = q_close(pairDisplayLists{k}(2)) - q_close(k);
        pairBottomMid = q_strip(k) + 1i*pairBottomRowY;
        focusBottom = pairBottomMid - 0.5*relBottom;
        neighbourBottom = pairBottomMid + 0.5*relBottom;
        allPanelCenters = [allPanelCenters; focusBottom; neighbourBottom; ...
            q_strip(k) + 1i*pairBottomRowY; q_strip(k) + 1i*verticalPlusBottomY; ...
            real(pairBottomMid) + 1i*bottomPairTitleY]; %#ok<AGROW>
        if showPairVelocityLabels
            allPanelCenters = [allPanelCenters; ...
                focusBottom - pairLabelXOffset - 1i*pairLabelYOffset; ...
                neighbourBottom + pairLabelXOffset - 1i*pairLabelYOffset]; %#ok<AGROW>
        end
    end
end

figure('Color','w','Name','mar26 visualise coarse fine - closest neighbours');
ax4 = axes;
setupFluidAxes(ax4,allPanelCenters,rad,fluidColor,1.7);
hold(ax4,'on');
annotationData.thirdFigureBounds = getBounds(allPanelCenters,rad,1.7);
for k = 1:numel(topDisplayCenters)
    psiText = ['\boldsymbol{\psi}^{(' num2str(k) ')}'];
    psiDisplayText = ['$\psi^{(' num2str(k) ')}$'];
    psiPos = [real(topDisplayCenters(k)), psiLabelY];
    annotationData.psiLabels(k).text = psiText;
    annotationData.psiLabels(k).position = psiPos;
    if drawAnnotations
        text(psiPos(1),psiPos(2),psiDisplayText, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);
    end
end

for k = 1:numel(topDisplayCenters)
    equalsPos = [real(topDisplayCenters(k)), equalsY];
    annotationData.equalsLabels(k).text = '=';
    annotationData.equalsLabels(k).position = equalsPos;
    if drawAnnotations
        text(equalsPos(1),equalsPos(2),'$=$', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Rotation',90,'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);
    end
end

numberLabelIdx = 1;
neighbourNumberLabelIdx = 1;
for k = 1:numel(topDisplayCenters)
    if k == 1 || k == numel(topDisplayCenters)
        drawDisk(ax4,topDisplayCenters(k),rad,bodyColor,'none',1.0,nCirclePts);
        drawBoundaryNodes(ax4,topDisplayCenters(k),rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
    end
    annotationData.circleNumberLabels(numberLabelIdx).text = sprintf('%d',k);
    annotationData.circleNumberLabels(numberLabelIdx).position = ...
        numberPosition(topDisplayCenters(k),annotationData.circleNumberLabels(numberLabelIdx).text);
    if drawAnnotations
        drawNumber(ax4,topDisplayCenters(k),annotationData.circleNumberLabels(numberLabelIdx).text,annotationTextColor);
    end
    numberLabelIdx = numberLabelIdx + 1;
end

for k = 1:numel(row2DisplayCenters)
    drawDisk(ax4,row2DisplayCenters(k),rad,bodyColor,'none',1.0,nCirclePts);
    drawBoundaryNodes(ax4,row2DisplayCenters(k),rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);

    phiText = ['\boldsymbol{\phi}^{(' num2str(k) ')}'];
    phiDisplayText = ['$\phi^{(' num2str(k) ')}$'];
    phiPos = [real(row2DisplayCenters(k)), phiLabelY];
    annotationData.phiLabels(k).text = phiText;
    annotationData.phiLabels(k).position = phiPos;
    if drawAnnotations
        text(phiPos(1),phiPos(2),phiDisplayText, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);
    end

    annotationData.circleNumberLabels(numberLabelIdx).text = sprintf('%d',k);
    annotationData.circleNumberLabels(numberLabelIdx).position = ...
        numberPosition(row2DisplayCenters(k),annotationData.circleNumberLabels(numberLabelIdx).text);
    if drawAnnotations
        drawNumber(ax4,row2DisplayCenters(k),annotationData.circleNumberLabels(numberLabelIdx).text,annotationTextColor);
    end
    numberLabelIdx = numberLabelIdx + 1;
end

verticalPlusIdx = 1;
pairFocusLabelIdx = 1;
pairNeighbourLabelIdx = 1;
pairTitleLabelIdx = 1;
for k = 1:numel(q_close)
    centerNow = q_strip(k) + 1i*topRowY;
    overallIdx = k + 1;
    neighboursNow = q_close(neighbourLists{k}) - q_close(k) + centerNow;

    for j = 1:numel(neighboursNow)
        drawDisk(ax4,neighboursNow(j),rad,bodyColor,blendWithWhite(bodyEdgeColor,0.45),0.9,nCirclePts);
        drawBoundaryNodes(ax4,neighboursNow(j),rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);
        neighbourOverallIdx = neighbourLists{k}(j) + 1;
        annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text = sprintf('%d',neighbourOverallIdx);
        annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).position = ...
            numberPosition(neighboursNow(j),annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text);
        if drawAnnotations
            drawNumber(ax4,neighboursNow(j), ...
                annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text, ...
                neighbourNumberColor);
        end
        neighbourNumberLabelIdx = neighbourNumberLabelIdx + 1;
    end

    drawDisk(ax4,centerNow,rad,bodyColor,bodyEdgeColor,1.1,nCirclePts);
    drawBoundaryNodes(ax4,centerNow,rad,fineNodeCount,fineNodeColor,fineNodeSize);
    annotationData.circleNumberLabels(numberLabelIdx).text = sprintf('%d',overallIdx);
    annotationData.circleNumberLabels(numberLabelIdx).position = ...
        numberPosition(centerNow,annotationData.circleNumberLabels(numberLabelIdx).text);
    if drawAnnotations
        drawNumber(ax4,centerNow,annotationData.circleNumberLabels(numberLabelIdx).text,annotationTextColor);
    end
    numberLabelIdx = numberLabelIdx + 1;

    verticalPlusTopPos = [real(q_strip(k)), verticalPlusTopY];
    annotationData.verticalPlusLabels(verticalPlusIdx).text = '+';
    annotationData.verticalPlusLabels(verticalPlusIdx).position = verticalPlusTopPos;
    verticalPlusIdx = verticalPlusIdx + 1;
    if drawAnnotations
        text(verticalPlusTopPos(1),verticalPlusTopPos(2),'$+$', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);
    end

    if ~isempty(pairDisplayLists{k})
        relTop = q_close(pairDisplayLists{k}(1)) - q_close(k);
        pairTopMid = q_strip(k) + 1i*pairTopRowY;
        focusTop = pairTopMid - 0.5*relTop;
        neighbourTop = pairTopMid + 0.5*relTop;
        neighbourOverallIdx = pairDisplayLists{k}(1) + 1;

        drawDisk(ax4,focusTop,rad,bodyColor,bodyEdgeColor,1.0,nCirclePts);
        drawBoundaryNodes(ax4,focusTop,rad,fineNodeCount,fineNodeColor,fineNodeSize);
        drawDisk(ax4,neighbourTop,rad,bodyColor,blendWithWhite(bodyEdgeColor,0.45),0.9,nCirclePts);
        drawBoundaryNodes(ax4,neighbourTop,rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);
        annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text = sprintf('%d',neighbourOverallIdx);
        annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).position = ...
            numberPosition(neighbourTop,annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text);
        if drawAnnotations
            drawNumber(ax4,neighbourTop,annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text, ...
                neighbourNumberColor);
        end
        neighbourNumberLabelIdx = neighbourNumberLabelIdx + 1;

        pairTopTitle = ['$\boldsymbol{\chi}^{(' num2str(overallIdx) ',' num2str(neighbourOverallIdx) ')}$'];
        pairTopTitleDisplay = ['$\chi^{(' num2str(overallIdx) ',' num2str(neighbourOverallIdx) ')}$'];
        pairTopTitlePos = [real(pairTopMid), topPairTitleY];
        annotationData.pairTitleLabels(pairTitleLabelIdx).text = pairTopTitle;
        annotationData.pairTitleLabels(pairTitleLabelIdx).position = pairTopTitlePos;
        if drawAnnotations
            text(pairTopTitlePos(1),pairTopTitlePos(2),pairTopTitleDisplay, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Interpreter','latex','FontSize',formulaFontSize, ...
                'Color',annotationTextColor,'Parent',ax4);
        end
        pairTitleLabelIdx = pairTitleLabelIdx + 1;

        if showPairVelocityLabels
            focusTopLabel = '$\mathbf{u} = \mathbf{0}$';
            focusTopPos = [real(focusTop) - pairLabelXOffset, imag(focusTop) - pairLabelYOffset];
            annotationData.pairFocusLabels(pairFocusLabelIdx).text = focusTopLabel;
            annotationData.pairFocusLabels(pairFocusLabelIdx).position = focusTopPos;
            if drawAnnotations
                text(focusTopPos(1),focusTopPos(2),focusTopLabel, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Interpreter','latex','FontSize',formulaFontSize, ...
                    'Color',pairAnnotationColor,'Parent',ax4);
            end
            pairFocusLabelIdx = pairFocusLabelIdx + 1;

            neighbourTopLabel = ['$\mathbf{u} = -\boldsymbol{\phi}^{(' num2str(overallIdx) ')}$'];
            neighbourTopPos = [real(neighbourTop) + pairLabelXOffset, imag(neighbourTop) - pairLabelYOffset];
            annotationData.pairNeighbourLabels(pairNeighbourLabelIdx).text = neighbourTopLabel;
            annotationData.pairNeighbourLabels(pairNeighbourLabelIdx).position = neighbourTopPos;
            if drawAnnotations
                text(neighbourTopPos(1),neighbourTopPos(2),neighbourTopLabel, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Interpreter','latex','FontSize',formulaFontSize, ...
                    'Color',pairAnnotationColor,'Parent',ax4);
            end
            pairNeighbourLabelIdx = pairNeighbourLabelIdx + 1;
        end

        annotationData.circleNumberLabels(numberLabelIdx).text = sprintf('%d',overallIdx);
        annotationData.circleNumberLabels(numberLabelIdx).position = ...
            numberPosition(focusTop,annotationData.circleNumberLabels(numberLabelIdx).text);
        if drawAnnotations
            drawNumber(ax4,focusTop,annotationData.circleNumberLabels(numberLabelIdx).text,annotationTextColor);
        end
        numberLabelIdx = numberLabelIdx + 1;
    end

    if numel(pairDisplayLists{k}) >= 2
        verticalPlusBottomPos = [real(q_strip(k)), verticalPlusBottomY];
        annotationData.verticalPlusLabels(verticalPlusIdx).text = '+';
        annotationData.verticalPlusLabels(verticalPlusIdx).position = verticalPlusBottomPos;
        verticalPlusIdx = verticalPlusIdx + 1;
        if drawAnnotations
            text(verticalPlusBottomPos(1),verticalPlusBottomPos(2),'$+$', ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Interpreter','latex','FontSize',formulaFontSize, ...
                'Color',annotationTextColor,'Parent',ax4);
        end

        relBottom = q_close(pairDisplayLists{k}(2)) - q_close(k);
        pairBottomMid = q_strip(k) + 1i*pairBottomRowY;
        focusBottom = pairBottomMid - 0.5*relBottom;
        neighbourBottom = pairBottomMid + 0.5*relBottom;
        neighbourOverallIdx = pairDisplayLists{k}(2) + 1;

        drawDisk(ax4,focusBottom,rad,bodyColor,bodyEdgeColor,1.0,nCirclePts);
        drawBoundaryNodes(ax4,focusBottom,rad,fineNodeCount,fineNodeColor,fineNodeSize);
        drawDisk(ax4,neighbourBottom,rad,bodyColor,blendWithWhite(bodyEdgeColor,0.45),0.9,nCirclePts);
        drawBoundaryNodes(ax4,neighbourBottom,rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);
        annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text = sprintf('%d',neighbourOverallIdx);
        annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).position = ...
            numberPosition(neighbourBottom,annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text);
        if drawAnnotations
            drawNumber(ax4,neighbourBottom,annotationData.neighbourCircleNumberLabels(neighbourNumberLabelIdx).text, ...
                neighbourNumberColor);
        end
        neighbourNumberLabelIdx = neighbourNumberLabelIdx + 1;

        pairBottomTitle = ['$\boldsymbol{\chi}^{(' num2str(overallIdx) ',' num2str(neighbourOverallIdx) ')}$'];
        pairBottomTitleDisplay = ['$\chi^{(' num2str(overallIdx) ',' num2str(neighbourOverallIdx) ')}$'];
        pairBottomTitlePos = [real(pairBottomMid), bottomPairTitleY];
        annotationData.pairTitleLabels(pairTitleLabelIdx).text = pairBottomTitle;
        annotationData.pairTitleLabels(pairTitleLabelIdx).position = pairBottomTitlePos;
        if drawAnnotations
            text(pairBottomTitlePos(1),pairBottomTitlePos(2),pairBottomTitleDisplay, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Interpreter','latex','FontSize',formulaFontSize, ...
                'Color',annotationTextColor,'Parent',ax4);
        end
        pairTitleLabelIdx = pairTitleLabelIdx + 1;

        if showPairVelocityLabels
            focusBottomLabel = '$\mathbf{u} = \mathbf{0}$';
            focusBottomPos = [real(focusBottom) - pairLabelXOffset, imag(focusBottom) - pairLabelYOffset];
            annotationData.pairFocusLabels(pairFocusLabelIdx).text = focusBottomLabel;
            annotationData.pairFocusLabels(pairFocusLabelIdx).position = focusBottomPos;
            if drawAnnotations
                text(focusBottomPos(1),focusBottomPos(2),focusBottomLabel, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Interpreter','latex','FontSize',formulaFontSize, ...
                    'Color',pairAnnotationColor,'Parent',ax4);
            end
            pairFocusLabelIdx = pairFocusLabelIdx + 1;

            neighbourBottomLabel = ['$\mathbf{u} = -\boldsymbol{\phi}^{(' num2str(overallIdx) ')}$'];
            neighbourBottomPos = [real(neighbourBottom) + pairLabelXOffset, imag(neighbourBottom) - pairLabelYOffset];
            annotationData.pairNeighbourLabels(pairNeighbourLabelIdx).text = neighbourBottomLabel;
            annotationData.pairNeighbourLabels(pairNeighbourLabelIdx).position = neighbourBottomPos;
            if drawAnnotations
                text(neighbourBottomPos(1),neighbourBottomPos(2),neighbourBottomLabel, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Interpreter','latex','FontSize',formulaFontSize, ...
                    'Color',pairAnnotationColor,'Parent',ax4);
            end
            pairNeighbourLabelIdx = pairNeighbourLabelIdx + 1;
        end

        annotationData.circleNumberLabels(numberLabelIdx).text = sprintf('%d',overallIdx);
        annotationData.circleNumberLabels(numberLabelIdx).position = ...
            numberPosition(focusBottom,annotationData.circleNumberLabels(numberLabelIdx).text);
        if drawAnnotations
            drawNumber(ax4,focusBottom,annotationData.circleNumberLabels(numberLabelIdx).text,annotationTextColor);
        end
        numberLabelIdx = numberLabelIdx + 1;
    end
end


if exportAnnotationData
    assignin('base',annotationVarName,annotationData);
end


function setupFluidAxes(ax,q,rad,fluidColor,padding)
bounds = getBounds(q,rad,padding);
hold(ax,'on');
patch(ax, ...
    [bounds(1), bounds(2), bounds(2), bounds(1)], ...
    [bounds(3), bounds(3), bounds(4), bounds(4)], ...
    fluidColor,'EdgeColor','none');
axis(ax,'equal');
xlim(ax,bounds(1:2));
ylim(ax,bounds(3:4));
axis(ax,'off');
set(ax,'Color','w');
end


function drawDisk(ax,center,rad,fillColor,edgeColor,lineWidth,nPts)
theta = linspace(0,2*pi,nPts+1);
theta(end) = [];
zc = center + rad*(cos(theta) + 1i*sin(theta));
patch(ax,real(zc),imag(zc),fillColor, ...
    'EdgeColor',edgeColor,'LineWidth',lineWidth);
end


function drawBoundaryNodes(ax,center,rad,nNodes,nodeColor,nodeSize)
theta = linspace(0,2*pi,nNodes+1);
theta(end) = [];
nodes = center + rad*(cos(theta) + 1i*sin(theta));
scatter(ax,real(nodes),imag(nodes),nodeSize,nodeColor,'filled', ...
    'MarkerEdgeColor',nodeColor,'LineWidth',0.7);
end


function drawCircleNumber(ax,center,labelText,labelColor,fontSize,fontName,baseOffset,sixOffset,sevenOffset)
centerPos = getCircleNumberPosition(center,labelText,baseOffset,sixOffset,sevenOffset);
text(centerPos(1),centerPos(2),labelText, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',fontSize,'FontName',fontName, ...
    'Color',labelColor,'Parent',ax);
end


function pos = getCircleNumberPosition(center,labelText,baseOffset,sixOffset,sevenOffset)
offset = baseOffset;
if strcmp(labelText,'6')
    offset = offset + sixOffset;
elseif strcmp(labelText,'7')
    offset = offset + sevenOffset;
end
pos = [real(center) + offset(1), imag(center) + offset(2)];
end


function bounds = getBounds(q,rad,padding)
xmin = min(real(q)) - rad - padding;
xmax = max(real(q)) + rad + padding;
ymin = min(imag(q)) - rad - padding;
ymax = max(imag(q)) + rad + padding;
bounds = [xmin, xmax, ymin, ymax];
end


function idx = closestNeighbours(q,k)
dist = abs(q - q(k));
dist(k) = inf;
minDist = min(dist);
idx = find(dist <= minDist + 1e-10*max(1,minDist));
end


function out = blendWithWhite(color,amount)
out = color + amount*(1 - color);
end
