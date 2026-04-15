function output = mar26_visualse_coarse_fine_new_rotation()
fprintf('=== Coarse/Fine Visualisation (Mar 27, 2026) - new rotation ===\n');

rad = 1;
domainPadding = 1.8;
nCirclePts = 240;

fluidColor = [0.82 0.91 0.98];
bodyColor = [1 1 1];
bodyEdgeColor = [0.72 0.77 0.82];
insetEdgeColor = blendWithWhite([0.30 0.42 0.55],0.25);

coarseNodeCount = 10;
coarseNodeColor = [0.10 0.34 0.82];
coarseNodeSize = 44;

fineNodeCount = 30;
fineNodeColor = [0.83 0.21 0.18];
fineNodeSize = 5;
fadedFineNodeColor = blendWithWhite(fineNodeColor,0.68);

formulaFontSize = 14;
numberFontSize = 12;
annotationTextColor = [0 0 0];
neighbourNumberColor = [0.55 0.55 0.55];
numberFontName = 'Courier';
% Nudge these offsets if individual numerals look slightly off-center.
numberBaseOffset = [0.00, 0.00];
numberSixOffset = [-0.05, 0.00];
numberSevenOffset = [-0.10, 0.00];
numberPosition = @(center,label) getCircleNumberPosition( ...
    center,label,numberBaseOffset,numberSixOffset,numberSevenOffset);

rotate_group = true;
showPairVelocityLabels = false;
% Keep row spacing uniform after rotation by leaving q_strip unshifted.
% Use topChiColumnGap and bottomChiColumnGap below to move the chi columns.
gapAfterEquals = 0.95;
gapAfterTopPlus = 1.8;
gapAfterBottomPlus = 1.9;
% Move the two chi-columns horizontally in the rotated layout.
topChiColumnGap = 23.92;
bottomChiColumnGap = 42.12;
pairTitleTopClearance = 0.45;
pairTitleBottomClearance = 0.60;
mainAxesPosition = [0.06, 0.36, 0.58, 0.56];
insetAxesPosition = [0.18, 0.07, 0.34, 0.22];
zoomPairLabels = [4 5];
figure5InsetGapInches = 0.42;
figure5InsetCenterOffsetInches = 0.00;
figure5InsetBottomMarginInches = 0.45;
zoomOutlineColor = blendWithWhite(coarseNodeColor,0.18);
zoomOutlineLineWidth = 1.15;
mainPanelPadX = 1.15;
mainPanelPadY = 1.80;
insetPadX = 1.20;
insetPadBottom = 1.20;
insetPadTop = 1.20;

closeGap = 0.15;
turnAngles = (pi/180)*[15, 70, -20, 55];
[q_close, q_far, q_fine] = buildFineGeometry(rad,closeGap,turnAngles,rotate_group);

panelSpacing = 6.7;
topRowY = 0;
psiLabelY = 4.2;
equalsY = -4.1;
coarseFocusRowY = -7.4 - gapAfterEquals;
row2LabelYOffset = 2.30;
phiLabelY = coarseFocusRowY + row2LabelYOffset;
pairTopRowY = -13.9 - gapAfterEquals - gapAfterTopPlus;
pairBottomRowY = -21.4 - gapAfterEquals - gapAfterTopPlus - gapAfterBottomPlus;
verticalPlusTopY = 0.44*(coarseFocusRowY + pairTopRowY);
verticalPlusBottomY = 0.47*(pairTopRowY + pairBottomRowY) + 0.2;
pairTitleYOffset = 1.10;
bottomPairTitleLift = 0.45;
isolatedPairPadX = 1.70;
isolatedPairPadY = 1.80;
isolatedPairTitleYOffset = 0.55;
isolatedFocusLabelOffset = [-1.95, -0.4];
isolatedNeighbourLabelOffset = [1.3, -1.20];
isolatedBodyLabelFontSize = 8;

q_strip = (0:numel(q_close)-1).' * panelSpacing;
topFarDisplays = [q_strip(1) - panelSpacing + 1i*topRowY; ...
    q_strip(end) + panelSpacing + 1i*topRowY];
topDisplayCenters = [topFarDisplays(1); q_strip + 1i*topRowY; topFarDisplays(2)];
row2DisplayCenters = topDisplayCenters + 1i*coarseFocusRowY;

neighbourLists = cell(numel(q_close),1);
pairDisplayLists = cell(numel(q_close),1);
for k = 1:numel(q_close)
    neighbourLists{k} = closestNeighbours(q_close,k);
    pairDisplayLists{k} = neighbourLists{k};
end

topPairTitleY = -inf;
bottomPairTitleY = -inf;
layoutRotate = @(z) -1i*conj(z);
allPanelCentersRot = [layoutRotate(topDisplayCenters(:)); ...
    layoutRotate(row2DisplayCenters(:)); ...
    layoutRotate(real(topDisplayCenters) + 1i*psiLabelY); ...
    layoutRotate(real(topDisplayCenters) + 1i*equalsY); ...
    layoutRotate(real(row2DisplayCenters) + 1i*phiLabelY)];
for k = 1:numel(q_close)
    centerNow = q_strip(k) + 1i*topRowY;
    centerDisplay = layoutRotate(centerNow);
    localOffsets = q_close(neighbourLists{k}) - q_close(k);
    localCenters = centerDisplay + localOffsets;
    allPanelCentersRot = [allPanelCentersRot; centerDisplay; localCenters(:); ... %#ok<AGROW>
        layoutRotate(q_strip(k) + 1i*verticalPlusTopY)];

    if ~isempty(pairDisplayLists{k})
        relTop = q_close(pairDisplayLists{k}(1)) - q_close(k);
        pairTopMid = layoutRotate(q_strip(k) + 1i*pairTopRowY);
        focusTop = pairTopMid - 0.5*relTop;
        neighbourTop = pairTopMid + 0.5*relTop;
        topPairTitleY = max(topPairTitleY, ...
            max(imag([focusTop, neighbourTop])) + rad + pairTitleYOffset);
        allPanelCentersRot = [allPanelCentersRot; focusTop; neighbourTop; pairTopMid]; %#ok<AGROW>
    end

    if numel(pairDisplayLists{k}) >= 2
        relBottom = q_close(pairDisplayLists{k}(2)) - q_close(k);
        pairBottomMid = layoutRotate(q_strip(k) + 1i*pairBottomRowY);
        focusBottom = pairBottomMid - 0.5*relBottom;
        neighbourBottom = pairBottomMid + 0.5*relBottom;
        bottomPairTitleY = max(bottomPairTitleY, ...
            max(imag([focusBottom, neighbourBottom])) + rad + pairTitleYOffset + bottomPairTitleLift);
        allPanelCentersRot = [allPanelCentersRot; focusBottom; neighbourBottom; ... %#ok<AGROW>
            pairBottomMid; layoutRotate(q_strip(k) + 1i*verticalPlusBottomY)];
    end
end
topPairTitleY = topPairTitleY + pairTitleTopClearance;
bottomPairTitleY = bottomPairTitleY + pairTitleBottomClearance;
topChiColumnX = verticalPlusTopY + topChiColumnGap;
bottomChiColumnX = verticalPlusBottomY + bottomChiColumnGap;

chiTitleCentersRot = [topChiColumnX - 1i*q_strip(:); ...
    bottomChiColumnX - 1i*q_strip(:)];
allPanelCentersRot = [allPanelCentersRot; chiTitleCentersRot];

fig4 = figure('Color','w','Name','mar26 visualise coarse fine new rotation - closest neighbours');
ax4 = axes;
setupFluidAxes(ax4,allPanelCentersRot,rad,fluidColor,1.7);
hold(ax4,'on');

numberLabelIdx = 1;
for k = 1:numel(topDisplayCenters)
    psiPos = layoutRotate(real(topDisplayCenters(k)) + 1i*psiLabelY);
    text(real(psiPos),imag(psiPos),['$\psi^{(' num2str(k) ')}$'], ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Interpreter','latex','FontSize',formulaFontSize, ...
        'Color',annotationTextColor,'Parent',ax4);

    equalsPos = layoutRotate(real(topDisplayCenters(k)) + 1i*equalsY);
    text(real(equalsPos),imag(equalsPos),'$=$', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Interpreter','latex','FontSize',formulaFontSize, ...
        'Color',annotationTextColor,'Parent',ax4);

    centerRot = layoutRotate(topDisplayCenters(k));
    if k == 1 || k == numel(topDisplayCenters)
        drawDisk(ax4,centerRot,rad,bodyColor,'none',1.0,nCirclePts);
        drawBoundaryNodes(ax4,centerRot,rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
    end
    drawCircleNumber(ax4,centerRot,sprintf('%d',k), ...
        annotationTextColor,numberFontSize,numberFontName, ...
        numberBaseOffset,numberSixOffset,numberSevenOffset);
    numberLabelIdx = numberLabelIdx + 1; %#ok<NASGU>
end

for k = 1:numel(row2DisplayCenters)
    centerRot = layoutRotate(row2DisplayCenters(k));
    drawDisk(ax4,centerRot,rad,bodyColor,'none',1.0,nCirclePts);
    drawBoundaryNodes(ax4,centerRot,rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);

    phiPos = layoutRotate(real(row2DisplayCenters(k)) + 1i*phiLabelY);
    text(real(phiPos),imag(phiPos),['$\phi^{(' num2str(k) ')}$'], ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Interpreter','latex','FontSize',formulaFontSize, ...
        'Color',annotationTextColor,'Parent',ax4);

    drawCircleNumber(ax4,centerRot,sprintf('%d',k), ...
        annotationTextColor,numberFontSize,numberFontName, ...
        numberBaseOffset,numberSixOffset,numberSevenOffset);
end

for k = 1:numel(q_close)
    centerNow = q_strip(k) + 1i*topRowY;
    centerRot = layoutRotate(centerNow);
    overallIdx = k + 1;
    neighbourOffsets = q_close(neighbourLists{k}) - q_close(k);

    for j = 1:numel(neighbourOffsets)
        neighbourRot = centerRot + neighbourOffsets(j);
        drawDisk(ax4,neighbourRot,rad,bodyColor,blendWithWhite(bodyEdgeColor,0.45),0.9,nCirclePts);
        drawBoundaryNodes(ax4,neighbourRot,rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);
        drawCircleNumber(ax4,neighbourRot, ...
            sprintf('%d',neighbourLists{k}(j) + 1),neighbourNumberColor, ...
            numberFontSize,numberFontName,numberBaseOffset,numberSixOffset,numberSevenOffset);
    end

    drawDisk(ax4,centerRot,rad,bodyColor,bodyEdgeColor,1.1,nCirclePts);
    drawBoundaryNodes(ax4,centerRot,rad,fineNodeCount,fineNodeColor,fineNodeSize);
    drawCircleNumber(ax4,centerRot,sprintf('%d',overallIdx), ...
        annotationTextColor,numberFontSize,numberFontName, ...
        numberBaseOffset,numberSixOffset,numberSevenOffset);

    plusTopPos = layoutRotate(q_strip(k) + 1i*verticalPlusTopY);
    text(real(plusTopPos),imag(plusTopPos),'$+$', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Interpreter','latex','FontSize',formulaFontSize, ...
        'Color',annotationTextColor,'Parent',ax4);

    if ~isempty(pairDisplayLists{k})
        relTop = q_close(pairDisplayLists{k}(1)) - q_close(k);
        pairTopMid = layoutRotate(q_strip(k) + 1i*pairTopRowY);
        focusTop = pairTopMid - 0.5*relTop;
        neighbourTop = pairTopMid + 0.5*relTop;
        neighbourOverallIdx = pairDisplayLists{k}(1) + 1;

        drawDisk(ax4,focusTop,rad,bodyColor,bodyEdgeColor,1.0,nCirclePts);
        drawBoundaryNodes(ax4,focusTop,rad,fineNodeCount,fineNodeColor,fineNodeSize);
        drawDisk(ax4,neighbourTop,rad,bodyColor,blendWithWhite(bodyEdgeColor,0.45),0.9,nCirclePts);
        drawBoundaryNodes(ax4,neighbourTop,rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);
        drawCircleNumber(ax4,focusTop,sprintf('%d',overallIdx), ...
            annotationTextColor,numberFontSize,numberFontName, ...
            numberBaseOffset,numberSixOffset,numberSevenOffset);
        drawCircleNumber(ax4,neighbourTop,sprintf('%d',neighbourOverallIdx), ...
            neighbourNumberColor,numberFontSize,numberFontName, ...
            numberBaseOffset,numberSixOffset,numberSevenOffset);

        pairTopTitlePos = topChiColumnX - 1i*q_strip(k);
        text(real(pairTopTitlePos),imag(pairTopTitlePos), ...
            ['$\chi^{(' num2str(overallIdx) ',' num2str(neighbourOverallIdx) ')}$'], ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);
    end

    if numel(pairDisplayLists{k}) >= 2
        plusBottomPos = layoutRotate(q_strip(k) + 1i*verticalPlusBottomY);
        text(real(plusBottomPos),imag(plusBottomPos),'$+$', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);

        relBottom = q_close(pairDisplayLists{k}(2)) - q_close(k);
        pairBottomMid = layoutRotate(q_strip(k) + 1i*pairBottomRowY);
        focusBottom = pairBottomMid - 0.5*relBottom;
        neighbourBottom = pairBottomMid + 0.5*relBottom;
        neighbourOverallIdx = pairDisplayLists{k}(2) + 1;

        drawDisk(ax4,focusBottom,rad,bodyColor,bodyEdgeColor,1.0,nCirclePts);
        drawBoundaryNodes(ax4,focusBottom,rad,fineNodeCount,fineNodeColor,fineNodeSize);
        drawDisk(ax4,neighbourBottom,rad,bodyColor,blendWithWhite(bodyEdgeColor,0.45),0.9,nCirclePts);
        drawBoundaryNodes(ax4,neighbourBottom,rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);
        drawCircleNumber(ax4,focusBottom,sprintf('%d',overallIdx), ...
            annotationTextColor,numberFontSize,numberFontName, ...
            numberBaseOffset,numberSixOffset,numberSevenOffset);
        drawCircleNumber(ax4,neighbourBottom,sprintf('%d',neighbourOverallIdx), ...
            neighbourNumberColor,numberFontSize,numberFontName, ...
            numberBaseOffset,numberSixOffset,numberSevenOffset);

        pairBottomTitlePos = bottomChiColumnX - 1i*q_strip(k);
        text(real(pairBottomTitlePos),imag(pairBottomTitlePos), ...
            ['$\chi^{(' num2str(overallIdx) ',' num2str(neighbourOverallIdx) ')}$'], ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Interpreter','latex','FontSize',formulaFontSize, ...
            'Color',annotationTextColor,'Parent',ax4);
    end
end

if showPairVelocityLabels
    warning('showPairVelocityLabels is not yet used in the rotated layout.');
end

geometryRotate = @(z) (-1i)*conj(z);
q_close_rot = geometryRotate(q_close);
q_far_rot = geometryRotate(q_far);
q_fine_rot = geometryRotate(q_fine);

closePairIdx = zoomPairLabels - 1;
if any(closePairIdx < 1) || any(closePairIdx > numel(q_close_rot))
    error('zoomPairLabels must refer to chain-body labels 2 through %d.',numel(q_close_rot) + 1);
end
pairCenters = q_close_rot(closePairIdx);

fig5 = figure('Color','w','Name','mar26 visualise coarse fine new rotation - coarse pair inset');
ax5 = axes('Parent',fig5,'Position',mainAxesPosition);
pairDisplayRotate = @(z) -z;
q_fine_pair_display = pairDisplayRotate(q_fine_rot);
pairCentersDisplay = pairDisplayRotate(pairCenters);

leftmostKeepCandidates = setdiff(1:numel(q_fine_pair_display),closePairIdx);
if isempty(leftmostKeepCandidates)
    leftmostMainRemoveIdx = [];
else
    [~, leftmostLocalIdx] = min(real(q_fine_pair_display(leftmostKeepCandidates)));
    leftmostMainRemoveIdx = leftmostKeepCandidates(leftmostLocalIdx);
end

mainKeepMask = true(numel(q_fine_pair_display),1);
if ~isempty(leftmostMainRemoveIdx)
    mainKeepMask(leftmostMainRemoveIdx) = false;
end

pairBounds = getBounds(pairCentersDisplay,rad,0.0);
insetLocalBounds = [pairBounds(1) - insetPadX, pairBounds(2) + insetPadX, ...
    pairBounds(3) - insetPadBottom, pairBounds(4) + insetPadTop];
mainPanelBounds = getBoundsXY(q_fine_pair_display(mainKeepMask),rad,mainPanelPadX,mainPanelPadY);

setupFluidAxesFromBounds(ax5,mainPanelBounds,fluidColor);
hold(ax5,'on');

for k = 1:numel(q_fine_pair_display)
    if ~mainKeepMask(k)
        continue;
    end
    drawDisk(ax5,q_fine_pair_display(k),rad,bodyColor,'none',1.0,nCirclePts);
    drawBoundaryNodes(ax5,q_fine_pair_display(k),rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
end
drawRectOutline(ax5,insetLocalBounds,zoomOutlineColor,zoomOutlineLineWidth);

ax5Inset = axes('Parent',fig5,'Position',insetAxesPosition);
hold(ax5Inset,'on');
drawFluidRect(ax5Inset,insetLocalBounds,fluidColor,'none',1.0);
axis(ax5Inset,'equal');
xlim(ax5Inset,insetLocalBounds(1:2));
ylim(ax5Inset,insetLocalBounds(3:4));
axis(ax5Inset,'off');
set(ax5Inset,'Color','w');

for k = 1:numel(q_fine_pair_display)
    centerNow = q_fine_pair_display(k);
    if ~diskIntersectsBounds(centerNow,rad,insetLocalBounds)
        continue;
    end

    isPairBody = ismember(k,closePairIdx);
    if isPairBody
        drawDisk(ax5Inset,centerNow,rad,bodyColor,'none',1.0,nCirclePts);
        drawBoundaryNodes(ax5Inset,centerNow,rad,fineNodeCount,fineNodeColor,fineNodeSize);
    else
        drawDisk(ax5Inset,centerNow,rad,fluidColor,'none',1.0,nCirclePts);
        drawBoundaryNodes(ax5Inset,centerNow,rad,coarseNodeCount,coarseNodeColor,coarseNodeSize);
    end
end

arrangeInsetAxesToMatchMainScale(fig5,ax5,mainPanelBounds,ax5Inset,insetLocalBounds, ...
    figure5InsetGapInches,figure5InsetCenterOffsetInches, ...
    figure5InsetBottomMarginInches);

output = struct();
output.q_close = q_close;
output.q_far = q_far;
output.q_fine = q_fine;
output.q_close_rot = q_close_rot;
output.q_far_rot = q_far_rot;
output.q_fine_rot = q_fine_rot;
output.rightmostPairIndices = closePairIdx;
output.figure4 = fig4;
output.figure4Axes = ax4;
output.figure5 = fig5;
output.figure5MainAxes = ax5;
output.figure5InsetAxes = ax5Inset;
output.figure5MainBounds = mainPanelBounds;
output.figure5InsetBounds = insetLocalBounds;
output.figure5InsetGapInches = figure5InsetGapInches;
output.figure5InsetCenterOffsetInches = figure5InsetCenterOffsetInches;
output.figure5InsetBottomMarginInches = figure5InsetBottomMarginInches;

pairOverallIdx = closePairIdx + 1;
%pairDistance = abs(pairCenters(2) - pairCenters(1))*0.7;
pairDistance = abs(pairCenters(2) - pairCenters(1))*0.93;
isolatedPairCenters = [(0.5i+0.1)*pairDistance; -(0.5i+0.1)*pairDistance];
%isolatedPairCenters = [(0.5i)*pairDistance; -(0.5i)*pairDistance];
isolatedPairBounds = getBoundsXY(isolatedPairCenters,rad,isolatedPairPadX,isolatedPairPadY);

fig6 = figure('Color','w','Name','mar26 visualise coarse fine new rotation - isolated pair');
ax6 = axes;
setupFluidAxesFromBounds(ax6,isolatedPairBounds,fluidColor);
hold(ax6,'on');

drawDisk(ax6,isolatedPairCenters(1),rad,bodyColor,'none',1.0,nCirclePts);
drawBoundaryNodes(ax6,isolatedPairCenters(1),rad,fineNodeCount,fineNodeColor,fineNodeSize);
drawDisk(ax6,isolatedPairCenters(2),rad,bodyColor,'none',0.9,nCirclePts);
drawBoundaryNodes(ax6,isolatedPairCenters(2),rad,fineNodeCount,fadedFineNodeColor,fineNodeSize);

isolatedPairTitlePos = [0, max(imag(isolatedPairCenters)) + rad + isolatedPairTitleYOffset];
text(isolatedPairTitlePos(1),isolatedPairTitlePos(2), ...
    ['$\mathbf{u} = \mathbf{0}$'], ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','latex','FontSize',formulaFontSize, ...
    'Color',annotationTextColor,'Parent',ax6);

text(real(isolatedPairCenters(1)),imag(isolatedPairCenters(1)),'$i$', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','latex','FontSize',isolatedBodyLabelFontSize, ...
    'Color',annotationTextColor,'Parent',ax6);
text(real(isolatedPairCenters(2)),imag(isolatedPairCenters(2)),'$j$', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'Interpreter','latex','FontSize',isolatedBodyLabelFontSize, ...
    'Color',annotationTextColor,'Parent',ax6);

focusDataPos = [real(isolatedPairCenters(1)) + isolatedFocusLabelOffset(1), ...
    imag(isolatedPairCenters(1)) + isolatedFocusLabelOffset(2)];
text(focusDataPos(1),focusDataPos(2),'$\chi^{(i,j)}$', ...
    'HorizontalAlignment','left','VerticalAlignment','middle', ...
    'Interpreter','latex','FontSize',formulaFontSize, ...
    'Color',annotationTextColor,'Parent',ax6);

neighbourDataPos = [real(isolatedPairCenters(2)) + isolatedNeighbourLabelOffset(1), ...
    imag(isolatedPairCenters(2)) + isolatedNeighbourLabelOffset(2)];
text(neighbourDataPos(1),neighbourDataPos(2), ...
    ['$\mathbf{u} = -\boldsymbol{\phi}^{(i)}$'], ...
    'HorizontalAlignment','right','VerticalAlignment','middle', ...
    'Interpreter','latex','FontSize',formulaFontSize, ...
    'Color',annotationTextColor,'Parent',ax6);

output.figure6 = fig6;
output.figure6Axes = ax6;
output.isolatedPairIndices = pairOverallIdx;
output.isolatedPairCenters = isolatedPairCenters;

end


function [q_close, q_far, q_fine] = buildFineGeometry(rad,closeGap,turnAngles,rotate_group)
dClose = rad*(2 + closeGap);
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
[q_close, ~] = sortBodiesTopToBottom(q_close);
[q_far, ~] = sortBodiesTopToBottom(q_far);
q_fine = [q_close; q_far];
end


function drawCircleNumber(ax,center,labelText,labelColor,fontSize,fontName,baseOffset,sixOffset,sevenOffset)
centerPos = getCircleNumberPosition(center,labelText,baseOffset,sixOffset,sevenOffset);
text(centerPos(1),centerPos(2),labelText, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',fontSize,'FontName',fontName, ...
    'Color',labelColor,'Parent',ax);
end


function drawCircleNumberRotated(ax,rotateFun,center,labelText,labelColor,fontSize,fontName,baseOffset,sixOffset,sevenOffset)
centerRot = rotateFun(center);
drawCircleNumber(ax,centerRot,labelText,labelColor,fontSize,fontName, ...
    baseOffset,sixOffset,sevenOffset);
end


function setupFluidAxes(ax,q,rad,fluidColor,padding)
bounds = getBounds(q,rad,padding);
setupFluidAxesFromBounds(ax,bounds,fluidColor);
end


function setupFluidAxesFromBounds(ax,bounds,fluidColor)
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


function drawFluidRect(ax,bounds,fillColor,edgeColor,lineWidth)
if ischar(edgeColor) || (isstring(edgeColor) && isscalar(edgeColor))
    edgeColorText = char(edgeColor);
else
    edgeColorText = '';
end

if strcmp(edgeColorText,'none')
    patch(ax, ...
        [bounds(1), bounds(2), bounds(2), bounds(1)], ...
        [bounds(3), bounds(3), bounds(4), bounds(4)], ...
        fillColor,'EdgeColor','none');
else
    patch(ax, ...
        [bounds(1), bounds(2), bounds(2), bounds(1)], ...
        [bounds(3), bounds(3), bounds(4), bounds(4)], ...
        fillColor,'EdgeColor',edgeColor,'LineWidth',lineWidth);
end
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
    set(fig,'Position',figPos);
    set(axMain,'Position',mainPos);
    insetY = bottomMarginInches;
end

set(axInset,'Position',[insetX, insetY, insetWidth, insetHeight]);
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


function restoreFigureAxesUnits(fig,axMain,axInset,figUnits,mainUnits,insetUnits)
set(fig,'Units',figUnits);
set(axMain,'Units',mainUnits);
set(axInset,'Units',insetUnits);
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


function bounds = getBoundsXY(q,rad,padX,padY)
xmin = min(real(q)) - rad - padX;
xmax = max(real(q)) + rad + padX;
ymin = min(imag(q)) - rad - padY;
ymax = max(imag(q)) + rad + padY;
bounds = [xmin, xmax, ymin, ymax];
end


function bounds = expandBounds(bounds,padding)
if isscalar(padding)
    padding = [padding, padding];
end
bounds = [bounds(1) - padding(1), bounds(2) + padding(1), ...
    bounds(3) - padding(2), bounds(4) + padding(2)];
end


function bounds = combineBounds(boundsA,boundsB,padding)
bounds = [min(boundsA(1),boundsB(1)), max(boundsA(2),boundsB(2)), ...
    min(boundsA(3),boundsB(3)), max(boundsA(4),boundsB(4))];
bounds = expandBounds(bounds,padding);
end


function drawRectOutline(ax,bounds,edgeColor,lineWidth)
patch(ax, ...
    [bounds(1), bounds(2), bounds(2), bounds(1)], ...
    [bounds(3), bounds(3), bounds(4), bounds(4)], ...
    'w','FaceColor','none','EdgeColor',edgeColor, ...
    'LineWidth',lineWidth,'LineStyle','-');
end


function drawInsetConnectors(fig,axMain,mainBounds,axInset,insetBounds,edgeColor,lineWidth)
drawnow;
mainTop = dataPointToFigure(axMain,[mainBounds(2), mainBounds(4)]);
mainBottom = dataPointToFigure(axMain,[mainBounds(2), mainBounds(3)]);
insetTop = dataPointToFigure(axInset,[insetBounds(1), insetBounds(4)]);
insetBottom = dataPointToFigure(axInset,[insetBounds(1), insetBounds(3)]);

axOverlay = axes('Parent',fig,'Position',[0 0 1 1], ...
    'Color','none','XLim',[0 1],'YLim',[0 1], ...
    'HitTest','off','HandleVisibility','off');
hold(axOverlay,'on');
line(axOverlay,[mainTop(1), insetTop(1)],[mainTop(2), insetTop(2)], ...
    'Color',edgeColor,'LineWidth',lineWidth,'Clipping','off');
line(axOverlay,[mainBottom(1), insetBottom(1)],[mainBottom(2), insetBottom(2)], ...
    'Color',edgeColor,'LineWidth',lineWidth,'Clipping','off');
axis(axOverlay,'off');
end


function figPoint = dataPointToFigure(ax,point)
axUnits = get(ax,'Units');
cleanupObj = onCleanup(@() set(ax,'Units',axUnits)); %#ok<NASGU>
set(ax,'Units','normalized');
axPos = get(ax,'Position');
xLim = xlim(ax);
yLim = ylim(ax);

u = (point(1) - xLim(1))/(xLim(2) - xLim(1));
v = (point(2) - yLim(1))/(yLim(2) - yLim(1));
figPoint = [axPos(1) + u*axPos(3), axPos(2) + v*axPos(4)];
end


function tf = diskIntersectsBounds(center,rad,bounds)
tf = real(center) + rad >= bounds(1) && real(center) - rad <= bounds(2) && ...
    imag(center) + rad >= bounds(3) && imag(center) - rad <= bounds(4);
end


function idx = rightmostClosePair(q_close_rot)
pairScores = -inf(numel(q_close_rot)-1,1);
for k = 1:numel(q_close_rot)-1
    pairScores(k) = real(mean(q_close_rot(k:k+1)));
end
[~, bestPair] = max(pairScores);
idx = bestPair:(bestPair+1);
end


function idx = closestNeighbours(q,k)
dist = abs(q - q(k));
dist(k) = inf;
minDist = min(dist);
idx = find(dist <= minDist + 1e-10*max(1,minDist));
end


function [qSorted, order] = sortBodiesTopToBottom(q)
sortData = [-imag(q(:)), real(q(:))];
[~, order] = sortrows(sortData,[1 2]);
qSorted = q(order);
end


function out = blendWithWhite(color,amount)
out = color + amount*(1 - color);
end
