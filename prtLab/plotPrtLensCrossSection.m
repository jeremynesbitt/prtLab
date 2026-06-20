function plotPrtLensCrossSection(T, options)
%PLOTPRTLENS cross-section plot for a prtLab optical system table.
%
%   plotPrtLensCrossSection(T) plots the y-z cross section of supported
%   plane and odd_asphere surfaces. Surface vertex positions are derived
%   from T.Thickness.

arguments
    T table
    options.NumPoints (1,1) double {mustBePositive} = 120
    options.LensColor (1,3) double = [0.60 0.80 1.00]
    options.FilterColor (1,3) double = [0.60 1.00 0.75]
    options.DefaultClearAperture (1,1) double {mustBePositive} = 1.0
    options.ShowLabels (1,1) logical = true
    options.ShowStop (1,1) logical = true
    options.StopZ (1,1) double = 0.20
    options.EntrancePupilRadius (1,1) double {mustBePositive} = 5.57/(2*2.8)
end

numSurfaces = height(T) - 1;
vertexZ = prtSurfaceVertexZ(T);

hold on;
axis equal;
box on;

surfaceIndex = 1;
while surfaceIndex <= numSurfaces
    if surfaceIndex < numSurfaces && isLensElementPair(T, surfaceIndex)
        plotElement(T, vertexZ, surfaceIndex, surfaceIndex+1, ...
            options.NumPoints, options.LensColor, options.DefaultClearAperture);
        surfaceIndex = surfaceIndex + 2;
    elseif surfaceIndex < numSurfaces && isPlanePair(T, surfaceIndex)
        plotPlanePair(T, vertexZ, surfaceIndex, surfaceIndex+1, ...
            options.NumPoints, options.FilterColor, options.DefaultClearAperture);
        surfaceIndex = surfaceIndex + 2;
    else
        plotSingleSurface(T, vertexZ, surfaceIndex, ...
            options.NumPoints, options.DefaultClearAperture);
        surfaceIndex = surfaceIndex + 1;
    end
end

if options.ShowStop
    stopExt = 0.25;
    plot([options.StopZ, options.StopZ], ...
        [options.EntrancePupilRadius, options.EntrancePupilRadius + stopExt], ...
        'k-', 'LineWidth', 3);
    plot([options.StopZ, options.StopZ], ...
        [-options.EntrancePupilRadius - stopExt, -options.EntrancePupilRadius], ...
        'k-', 'LineWidth', 3);
end

zMin = min(vertexZ(1:numSurfaces)) - 0.5;
zMax = max(vertexZ(1:numSurfaces)) + 0.7;
plot([zMin, zMax], [0, 0], 'k:', 'LineWidth', 0.8);

if options.ShowLabels
    for surfaceIndex = 1:numSurfaces
        clearAperture = prtClearAperture(T(surfaceIndex,:), options.DefaultClearAperture);
        zLabel = vertexZ(surfaceIndex) + prtSurfaceSag(clearAperture*0.85, 0, T(surfaceIndex,:));
        text(zLabel, clearAperture + 0.12, num2str(surfaceIndex), ...
            'FontSize', 9, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

xlabel('z (mm)');
ylabel('y (mm)');
xlim([zMin, zMax]);

end

function plotElement(T, vertexZ, indexA, indexB, numPoints, color, defaultClearAperture)
clearAperture = min( ...
    prtClearAperture(T(indexA,:), defaultClearAperture), ...
    prtClearAperture(T(indexB,:), defaultClearAperture));
y = linspace(-clearAperture, clearAperture, 2*numPoints - 1);
zA = arrayfun(@(yy) vertexZ(indexA) + prtSurfaceSag(0, yy, T(indexA,:)), y);
zB = arrayfun(@(yy) vertexZ(indexB) + prtSurfaceSag(0, yy, T(indexB,:)), y);
fill([zA, fliplr(zB)], [y, fliplr(y)], color, ...
    'EdgeColor', 'k', 'LineWidth', 1.2);
end

function plotPlanePair(T, vertexZ, indexA, indexB, numPoints, color, defaultClearAperture)
clearAperture = min( ...
    prtClearAperture(T(indexA,:), defaultClearAperture), ...
    prtClearAperture(T(indexB,:), defaultClearAperture));
y = linspace(-clearAperture, clearAperture, numPoints);
fill([vertexZ(indexA)*ones(1,numPoints), vertexZ(indexB)*ones(1,numPoints)], ...
    [y, fliplr(y)], color, 'EdgeColor', 'k', 'LineWidth', 1.2);
end

function plotSingleSurface(T, vertexZ, surfaceIndex, numPoints, defaultClearAperture)
clearAperture = prtClearAperture(T(surfaceIndex,:), defaultClearAperture);
y = linspace(-clearAperture, clearAperture, numPoints);
z = arrayfun(@(yy) vertexZ(surfaceIndex) + prtSurfaceSag(0, yy, T(surfaceIndex,:)), y);
plot(z, y, 'k-', 'LineWidth', 1.2);
end

function tf = isLensElementPair(T, surfaceIndex)
if T.MaterialType(surfaceIndex+1) ~= "isotropic"
    tf = false;
    return;
end

nNext = T.IndexData{surfaceIndex+1}.n;
tf = isreal(nNext) && nNext ~= 1;
end

function tf = isPlanePair(T, surfaceIndex)
tf = string(T.SurfaceType(surfaceIndex)) == "plane" && ...
    string(T.SurfaceType(surfaceIndex+1)) == "plane";
end
