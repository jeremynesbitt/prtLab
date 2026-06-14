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
end

hold on;
outputs = normalizeRayOutputs(rayOutputs);
colors = normalizeColors(options.Color, numel(outputs));

for outputIndex = 1:numel(outputs)
    rayOutput = outputs{outputIndex};
    finalRayIds = rayOutput.finalRayIds;
    if isempty(finalRayIds) && options.PlotIncomplete
        plotIncompleteRay(rayOutput, colors(outputIndex,:), options.LineWidth, ...
            options.PreExtend, options.IncompleteLineStyle);
    else
        for finalIndex = 1:numel(finalRayIds)
            plotOneFinalRay(rayOutput, finalRayIds(finalIndex), ...
            colors(outputIndex,:), options.LineWidth, ...
            options.PreExtend, options.PostExtend);
        end
    end
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

function plotOneFinalRay(rayOutput, finalRayId, color, lineWidth, preExtend, postExtend)
history = rayOutput.rays(finalRayId).history;
positions = zeros(numel(history), 3);
for historyIndex = 1:numel(history)
    positions(historyIndex,:) = rayOutput.rays(history(historyIndex)).position(:).';
end

initialRay = rayOutput.rays(history(1));
initialPoint = positions(1,:);
initialDirection = initialRay.S(:).';
prePoint = initialPoint - preExtend*initialDirection;

plot([prePoint(3), initialPoint(3)], [prePoint(2), initialPoint(2)], ...
    '-', 'Color', color, 'LineWidth', lineWidth);
plot(positions(:,3), positions(:,2), ...
    '-', 'Color', color, 'LineWidth', lineWidth);

finalRay = rayOutput.rays(finalRayId);
finalPoint = positions(end,:);
finalDirection = finalRay.S(:).';
postPoint = finalPoint + postExtend*finalDirection;
plot([finalPoint(3), postPoint(3)], [finalPoint(2), postPoint(2)], ...
    '-', 'Color', color, 'LineWidth', lineWidth);
end

function plotIncompleteRay(rayOutput, color, lineWidth, preExtend, lineStyle)
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
end
