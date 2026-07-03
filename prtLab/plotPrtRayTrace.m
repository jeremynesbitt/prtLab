function plotPrtRayTrace(rayOutputs, options)
%PLOTPRTRAYTRACE Overlay polarizationRayTrace ray paths in y-z.
%
%   plotPrtRayTrace(rayOutput)
%   plotPrtRayTrace({rayOutput1, rayOutput2, ...})
%
%   For branched traces, each final ray history is plotted.

arguments
    rayOutputs
    options.Color = [0.85 0.10 0.10]
    options.LineWidth (1,1) double {mustBePositive} = 1.5
    options.PreExtend (1,1) double {mustBeNonnegative} = 0.3
    options.PostExtend (1,1) double {mustBeNonnegative} = 0.5
    options.PlotIncomplete (1,1) logical = true
    options.IncompleteLineStyle (1,1) string = "--"
    options.AutoExpandAxes (1,1) logical = true
    options.AxisPaddingFraction (1,1) double {mustBeNonnegative} = 0.04
    options.ColorByMode (1,1) logical = true
    options.ModeColors struct = struct( ...
        'input', [0.85 0.10 0.10], ...
        'transmitted', [0.85 0.10 0.10], ...
        'ordinary', [0.00 0.30 0.90], ...
        'extraordinary', [0.00 0.55 0.20], ...
        'reflected', [0.95 0.45 0.05])
end

hold on;
outputs = normalizeRayOutputs(rayOutputs);
colors = normalizeColors(options.Color, numel(outputs));
plottedPoints = zeros(0, 2);

for outputIndex = 1:numel(outputs)
    rayOutput = outputs{outputIndex};
    finalRayIds = rayOutput.finalRayIds;
    if isempty(finalRayIds) && options.PlotIncomplete
        points = plotIncompleteRay(rayOutput, colors(outputIndex,:), options.LineWidth, ...
            options.PreExtend, options.IncompleteLineStyle);
        plottedPoints = [plottedPoints; points]; %#ok<AGROW>
    else
        for finalIndex = 1:numel(finalRayIds)
            points = plotOneFinalRay(rayOutput, finalRayIds(finalIndex), ...
            colors(outputIndex,:), options.LineWidth, ...
            options.PreExtend, options.PostExtend, ...
            options.ColorByMode, options.ModeColors);
            plottedPoints = [plottedPoints; points]; %#ok<AGROW>
        end
    end
end

if options.AutoExpandAxes
    expandAxesToPoints(plottedPoints, options.AxisPaddingFraction);
end

end

function outputs = normalizeRayOutputs(rayOutputs)
if iscell(rayOutputs)
    outputs = rayOutputs(:).';
elseif numel(rayOutputs) > 1
    outputs = num2cell(rayOutputs);
else
    outputs = {rayOutputs};
end
end

function colors = normalizeColors(colorInput, numOutputs)
if ischar(colorInput) || isstring(colorInput)
    colors = repmat(colorInput, numOutputs, 1);
    return;
end

if size(colorInput, 1) == 1
    colors = repmat(colorInput, numOutputs, 1);
else
    colors = colorInput;
end
end

function points = plotOneFinalRay(rayOutput, finalRayId, color, lineWidth, preExtend, postExtend, colorByMode, modeColors)
history = rayOutput.rays(finalRayId).history;
positions = zeros(numel(history), 3);
for historyIndex = 1:numel(history)
    positions(historyIndex,:) = rayOutput.rays(history(historyIndex)).position(:).';
end

initialRay = rayOutput.rays(history(1));
initialPoint = positions(1,:);
initialDirection = initialRay.S(:).';
prePoint = initialPoint - preExtend*initialDirection;

preColor = rayColor(initialRay, color, colorByMode, modeColors);
plot([prePoint(3), initialPoint(3)], [prePoint(2), initialPoint(2)], ...
    '-', 'Color', preColor, 'LineWidth', lineWidth);

for segmentIndex = 2:numel(history)
    segmentRay = rayOutput.rays(history(segmentIndex-1));
    segmentColor = rayColor(segmentRay, color, colorByMode, modeColors);
    plot(positions(segmentIndex-1:segmentIndex,3), positions(segmentIndex-1:segmentIndex,2), ...
        '-', 'Color', segmentColor, 'LineWidth', lineWidth);
end

finalRay = rayOutput.rays(finalRayId);
finalPoint = positions(end,:);
finalDirection = finalRay.S(:).';
postPoint = finalPoint + postExtend*finalDirection;
postColor = rayColor(finalRay, color, colorByMode, modeColors);
plot([finalPoint(3), postPoint(3)], [finalPoint(2), postPoint(2)], ...
    '-', 'Color', postColor, 'LineWidth', lineWidth);

points3 = [prePoint; positions; postPoint];
points = points3(:,[3, 2]);
end

function color = rayColor(ray, fallbackColor, colorByMode, modeColors)
color = fallbackColor;
if ~colorByMode || ~isfield(ray, 'mode')
    return;
end

modeName = char(ray.mode);
if isfield(modeColors, modeName)
    color = modeColors.(modeName);
end
end

function points = plotIncompleteRay(rayOutput, color, lineWidth, preExtend, lineStyle)
points = zeros(0, 2);
if isempty(rayOutput.rays)
    return;
end

positions = zeros(numel(rayOutput.interactions)+1, 3);
positions(1,:) = rayOutput.rays(1).position(:).';
for interactionIndex = 1:numel(rayOutput.interactions)
    positions(interactionIndex+1,:) = ...
        rayOutput.interactions(interactionIndex).position(:).';
end

initialRay = rayOutput.rays(1);
initialPoint = positions(1,:);
initialDirection = initialRay.S(:).';
prePoint = initialPoint - preExtend*initialDirection;

plot([prePoint(3), initialPoint(3)], [prePoint(2), initialPoint(2)], ...
    lineStyle, 'Color', color, 'LineWidth', lineWidth);
plot(positions(:,3), positions(:,2), ...
    lineStyle, 'Color', color, 'LineWidth', lineWidth);

points3 = [prePoint; positions];
points = points3(:,[3, 2]);
end

function expandAxesToPoints(points, paddingFraction)
if isempty(points)
    return;
end

ax = gca;
currentXLim = xlim(ax);
currentYLim = ylim(ax);

newXLim = [min([currentXLim(1); points(:,1)]), ...
           max([currentXLim(2); points(:,1)])];
newYLim = [min([currentYLim(1); points(:,2)]), ...
           max([currentYLim(2); points(:,2)])];

xPadding = paddingFraction*max(eps, diff(newXLim));
yPadding = paddingFraction*max(eps, diff(newYLim));

xlim(ax, newXLim + [-xPadding, xPadding]);
ylim(ax, newYLim + [-yPadding, yPadding]);
end
