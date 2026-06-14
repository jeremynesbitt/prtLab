function [T, clearApertures] = updateClearAperturesFromRayTrace(T, rayOutputs, options)
%UPDATECLEARAPERTURESFROMRAYTRACE Set plot clear apertures from traced rays.
%
%   T = updateClearAperturesFromRayTrace(T, rayOutput)
%   T = updateClearAperturesFromRayTrace(T, {rayOutput1, rayOutput2, ...})
%   [T, clearApertures] = updateClearAperturesFromRayTrace(...)
%
%   The clear aperture for each traced surface is set from the maximum
%   radial intercept sqrt(x^2+y^2) across all supplied ray outputs. This is
%   intended for plotting only; polarizationRayTrace does not clip rays by
%   clear aperture.

arguments
    T table
    rayOutputs
    options.Margin (1,1) double {mustBeNonnegative} = 1.0
    options.Minimum (1,1) double {mustBeNonnegative} = 0
    options.Surfaces double = 1:height(T)-1
end

numSurfaces = height(T) - 1;
clearApertures = zeros(numSurfaces, 1);
outputs = normalizeRayOutputs(rayOutputs);

for outputIndex = 1:numel(outputs)
    rayOutput = outputs{outputIndex};
    for interactionIndex = 1:numel(rayOutput.interactions)
        surfaceIndex = rayOutput.interactions(interactionIndex).surfaceIndex;
        if surfaceIndex < 1 || surfaceIndex > numSurfaces
            continue;
        end
        hit = rayOutput.interactions(interactionIndex).position;
        radialHeight = hypot(hit(1), hit(2));
        clearApertures(surfaceIndex) = max(clearApertures(surfaceIndex), radialHeight);
    end
end

clearApertures = max(options.Margin * clearApertures, options.Minimum);
surfacesToUpdate = intersect(options.Surfaces(:).', 1:numSurfaces);

for surfaceIndex = surfacesToUpdate
    surfaceData = T.SurfaceData{surfaceIndex};
    surfaceData.clearAperture = clearApertures(surfaceIndex);
    T.SurfaceData{surfaceIndex} = surfaceData;
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
