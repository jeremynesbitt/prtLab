function h = plotPrtSystem3D(T, rayOutputs, options)
%PLOTPRTSYSTEM3D Plot a prtLab system, ray paths, and polarization glyphs.
%
%   h = plotPrtSystem3D(T, rayOutput) renders supported surfaces, traced ray
%   branches, and small polarization ellipses from each ray branch field.
%
%   This is a branch-level visualization. For split o/e rays, the ordinary
%   and extraordinary branch fields are shown separately.

arguments
    T table
    rayOutputs
    options.ShowSystem (1,1) logical = true
    options.ShowRays (1,1) logical = true
    options.ShowPolarization (1,1) logical = true
    options.PolarizationAt (1,1) string {mustBeMember(options.PolarizationAt, ["vertices", "interfaces", "final", "all"])} = "interfaces"
    options.FieldName (1,1) string {mustBeMember(options.FieldName, ["fieldE", "modeE"])} = "fieldE"
    options.PolarizationScale = "auto"
    options.PolarizationScaleFraction (1,1) double {mustBePositive} = 0.10
    options.NormalizePolarization (1,1) logical = true
    options.SamplesPerEllipse (1,1) double {mustBeInteger, mustBePositive} = 101
    options.SurfaceSamples (1,1) double {mustBeInteger, mustBePositive} = 36
    options.DefaultClearAperture (1,1) double {mustBePositive} = 1.0
    options.SurfaceAlpha (1,1) double {mustBeGreaterThanOrEqual(options.SurfaceAlpha,0), mustBeLessThanOrEqual(options.SurfaceAlpha,1)} = 0.18
    options.EdgeAlpha (1,1) double {mustBeGreaterThanOrEqual(options.EdgeAlpha,0), mustBeLessThanOrEqual(options.EdgeAlpha,1)} = 0.06
    options.LineWidth (1,1) double {mustBePositive} = 1.5
    options.GlyphLineWidth (1,1) double {mustBePositive} = 1.4
    options.PreExtend (1,1) double {mustBeNonnegative} = 0.0
    options.PostExtend (1,1) double {mustBeNonnegative} = 0.0
    options.AxisPaddingFraction (1,1) double {mustBeNonnegative} = 0.16
    options.EqualAxes (1,1) logical = true
    options.ShowLegend (1,1) logical = false
    options.ColorByMode (1,1) logical = true
    options.ModeColors struct = struct( ...
        'input', [0.12 0.12 0.12], ...
        'transmitted', [0.12 0.12 0.12], ...
        'ordinary', [0.00 0.30 0.90], ...
        'extraordinary', [0.00 0.55 0.20], ...
        'reflected', [0.95 0.45 0.05])
    options.PolarizationColor (1,3) double = [0.85 0.10 0.10]
    options.View double = [35 20]
end

hold on;
if options.EqualAxes
    axis equal;
end
grid on;
box on;

outputs = normalizeRayOutputs(rayOutputs);
allPoints = zeros(0,3);
options.PolarizationScale = resolvePolarizationScale(T, outputs, options);
h = initializeHandles();

if options.ShowSystem
    [systemPoints, systemHandles] = plotSystemSurfaces(T, options);
    allPoints = [allPoints; systemPoints];
    h.surfaces = [h.surfaces; systemHandles(:)];
end

if options.ShowRays || options.ShowPolarization
    for outputIndex = 1:numel(outputs)
        rayOutput = outputs{outputIndex};
        if isempty(rayOutput.finalRayIds)
            continue;
        end
        for finalIndex = 1:numel(rayOutput.finalRayIds)
            finalRayId = rayOutput.finalRayIds(finalIndex);
            history = rayOutput.rays(finalRayId).history;
            [rayPoints, rayHandles] = plotRayHistory(rayOutput, history, options);
            allPoints = [allPoints; rayPoints];
            h.rays = [h.rays; rayHandles(:)];
            if options.ShowPolarization
                glyphHandles = plotPolarizationForHistory(rayOutput, history, options);
                h.glyphs = [h.glyphs; glyphHandles(:)];
            end
        end
    end
end

units = prtLengthUnits(T);
xlabel("x (" + units + ")");
ylabel("y (" + units + ")");
zlabel("z (" + units + ")");
view(options.View);

if ~isempty(allPoints)
    padAxes3(allPoints, options.AxisPaddingFraction);
end

if options.ShowLegend
    addModeLegend(options);
end

end

function h = initializeHandles()
h = struct( ...
    'surfaces', gobjects(0), ...
    'rays', gobjects(0), ...
    'glyphs', gobjects(0));
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

function [points, handles] = plotSystemSurfaces(T, options)
numSurfaces = height(T) - 1;
vertexZ = prtSurfaceVertexZ(T);
points = zeros(0,3);
handles = gobjects(0);

surfaceIndex = 1;
while surfaceIndex <= numSurfaces
    if surfaceIndex < numSurfaces && isPlanePair(T, surfaceIndex) && isMaterialRegion(T, surfaceIndex)
        [newPoints, newHandles] = plotPlanePairVolume(T, vertexZ, surfaceIndex, surfaceIndex+1, options);
        points = [points; newPoints];
        handles = [handles; newHandles(:)];
        surfaceIndex = surfaceIndex + 2;
        continue;
    end

    clearAperture = prtClearAperture(T(surfaceIndex,:), options.DefaultClearAperture);
    [xx, yy, zz] = surfaceMesh(T, vertexZ, surfaceIndex, clearAperture, options.SurfaceSamples);
    newHandle = plotSurfaceMesh(xx, yy, zz, surfaceColor(T, surfaceIndex), options.SurfaceAlpha, options.EdgeAlpha);
    points = [points; [xx(:), yy(:), zz(:)]];
    handles = [handles; newHandle];
    surfaceIndex = surfaceIndex + 1;
end
end

function tf = isPlanePair(T, surfaceIndex)
tf = string(T.SurfaceType(surfaceIndex)) == "plane" && ...
    string(T.SurfaceType(surfaceIndex+1)) == "plane";
end

function tf = isMaterialRegion(T, surfaceIndex)
materialIndex = surfaceIndex + 1;
tf = materialIndex <= height(T) && T.MaterialType(materialIndex) ~= "isotropic";
if tf
    return;
end

if materialIndex <= height(T) && T.MaterialType(materialIndex) == "isotropic"
    indexData = T.IndexData{materialIndex};
    tf = isfield(indexData, 'n') && indexData.n ~= 1;
end
end

function [points, handles] = plotPlanePairVolume(T, vertexZ, indexA, indexB, options)
clearAperture = min( ...
    prtClearAperture(T(indexA,:), options.DefaultClearAperture), ...
    prtClearAperture(T(indexB,:), options.DefaultClearAperture));
color = surfaceColor(T, indexA);
[xA, yA, zA] = surfaceMesh(T, vertexZ, indexA, clearAperture, options.SurfaceSamples);
[xB, yB, zB] = surfaceMesh(T, vertexZ, indexB, clearAperture, options.SurfaceSamples);

hA = plotSurfaceMesh(xA, yA, zA, color, options.SurfaceAlpha, options.EdgeAlpha);
hB = plotSurfaceMesh(xB, yB, zB, color, options.SurfaceAlpha, options.EdgeAlpha);

phi = linspace(0, 2*pi, options.SurfaceSamples);
xSide = clearAperture * [cos(phi); cos(phi)];
ySide = clearAperture * [sin(phi); sin(phi)];
zSide = [vertexZ(indexA) * ones(size(phi)); vertexZ(indexB) * ones(size(phi))];
hSide = surf(xSide, ySide, zSide, ...
    'FaceColor', color, ...
    'FaceAlpha', min(0.35, 1.7*options.SurfaceAlpha), ...
    'EdgeAlpha', options.EdgeAlpha, ...
    'EdgeColor', [0 0 0]);

points = [xA(:), yA(:), zA(:); xB(:), yB(:), zB(:); xSide(:), ySide(:), zSide(:)];
handles = [hA; hB; hSide];
end

function [xx, yy, zz] = surfaceMesh(T, vertexZ, surfaceIndex, clearAperture, numSamples)
r = linspace(0, clearAperture, numSamples);
phi = linspace(0, 2*pi, numSamples);
[rr, pp] = meshgrid(r, phi);
xx = rr .* cos(pp);
yy = rr .* sin(pp);
zz = zeros(size(xx));
for ii = 1:numel(xx)
    zz(ii) = vertexZ(surfaceIndex) + prtSurfaceSag(xx(ii), yy(ii), T(surfaceIndex,:));
end
end

function hSurface = plotSurfaceMesh(xx, yy, zz, color, alphaValue, edgeAlpha)
hSurface = surf(xx, yy, zz, ...
    'FaceColor', color, ...
    'FaceAlpha', alphaValue, ...
    'EdgeAlpha', edgeAlpha, ...
    'EdgeColor', [0 0 0]);
end

function color = surfaceColor(T, surfaceIndex)
if surfaceIndex + 1 <= height(T) && T.MaterialType(surfaceIndex+1) == "uniaxial"
    color = [0.60 1.00 0.75];
elseif surfaceIndex + 1 <= height(T) && T.MaterialType(surfaceIndex+1) == "isotropic"
    indexData = T.IndexData{surfaceIndex+1};
    if isfield(indexData, 'n') && isreal(indexData.n) && indexData.n ~= 1
        color = [0.60 0.80 1.00];
    else
        color = [0.85 0.85 0.85];
    end
else
    color = [0.85 0.85 0.85];
end
end

function [points, handles] = plotRayHistory(rayOutput, history, options)
positions = historyPositions(rayOutput, history);
points = positions;
handles = gobjects(0);

if options.PreExtend > 0
    initialRay = rayOutput.rays(history(1));
    prePoint = positions(1,:) - options.PreExtend * initialRay.S(:).';
    hLine = plot3([prePoint(1), positions(1,1)], [prePoint(2), positions(1,2)], [prePoint(3), positions(1,3)], ...
        '-', 'Color', rayColor(initialRay, options), 'LineWidth', options.LineWidth);
    points = [prePoint; points];
    handles = [handles; hLine];
end

for segmentIndex = 2:numel(history)
    segmentRay = rayOutput.rays(history(segmentIndex-1));
    hLine = plot3(positions(segmentIndex-1:segmentIndex,1), ...
        positions(segmentIndex-1:segmentIndex,2), ...
        positions(segmentIndex-1:segmentIndex,3), ...
        '-', 'Color', rayColor(segmentRay, options), 'LineWidth', options.LineWidth);
    handles = [handles; hLine];
end

if options.PostExtend > 0
    finalRay = rayOutput.rays(history(end));
    postPoint = positions(end,:) + options.PostExtend * finalRay.S(:).';
    hLine = plot3([positions(end,1), postPoint(1)], [positions(end,2), postPoint(2)], [positions(end,3), postPoint(3)], ...
        '-', 'Color', rayColor(finalRay, options), 'LineWidth', options.LineWidth);
    points = [points; postPoint];
    handles = [handles; hLine];
end
end

function positions = historyPositions(rayOutput, history)
positions = zeros(numel(history), 3);
for historyIndex = 1:numel(history)
    positions(historyIndex,:) = rayOutput.rays(history(historyIndex)).position(:).';
end
end

function handles = plotPolarizationForHistory(rayOutput, history, options)
switch options.PolarizationAt
    case "interfaces"
        rayIds = unique(history, 'stable');
    case "vertices"
        rayIds = history;
    case "final"
        rayIds = history(end);
    case "all"
        rayIds = history;
end

handles = gobjects(0);
plottedKeys = strings(0);
for rayId = rayIds(:).'
    ray = rayOutput.rays(rayId);
    key = sprintf('%.12g_%.12g_%.12g_%s', ray.position(1), ray.position(2), ray.position(3), ray.mode);
    if options.PolarizationAt == "interfaces" && any(plottedKeys == key)
        continue;
    end
    plottedKeys(end+1) = key; %#ok<AGROW>
    field = ray.(options.FieldName);
    color = options.PolarizationColor;
    glyphHandles = plotPolarizationGlyph3D(ray.position, field, color, options);
    handles = [handles; glyphHandles(:)];
end
end

function handles = plotPolarizationGlyph3D(center, field, color, options)
handles = gobjects(0);
if norm(field) < eps
    return;
end

if options.NormalizePolarization
    field = field / norm(field);
end
field = options.PolarizationScale * field;

theta = linspace(0, 2*pi, options.SamplesPerEllipse);
ellipse = real(exp(-1i*theta) .* field(:));
points = center(:) + ellipse;
hEllipse = plot3(points(1,:), points(2,:), points(3,:), ...
    '-', 'Color', color, 'LineWidth', options.GlyphLineWidth);
handles = [handles; hEllipse];

theta0 = -pi/6;
dtheta = 1e-2;
p0 = center(:) + real(exp(-1i*theta0) * field(:));
p1 = center(:) + real(exp(-1i*(theta0 + dtheta)) * field(:));
tangent = p1 - p0;
if norm(tangent) < eps
    return;
end
tangent = tangent / norm(tangent);
hArrow = quiver3(p0(1), p0(2), p0(3), ...
    0.35*options.PolarizationScale*tangent(1), ...
    0.35*options.PolarizationScale*tangent(2), ...
    0.35*options.PolarizationScale*tangent(3), ...
    0, 'Color', color, 'LineWidth', 1.0, 'MaxHeadSize', 1.5);
handles = [handles; hArrow];
end

function color = rayColor(ray, options)
color = [0.85 0.10 0.10];
if ~options.ColorByMode || ~isfield(ray, 'mode')
    return;
end

modeName = char(ray.mode);
if isfield(options.ModeColors, modeName)
    color = options.ModeColors.(modeName);
end
end

function scale = resolvePolarizationScale(T, outputs, options)
scale = options.PolarizationScale;
if isnumeric(scale)
    if ~isscalar(scale) || scale <= 0
        error('plotPrtSystem3D:InvalidPolarizationScale', ...
            'PolarizationScale must be "auto" or a positive scalar.');
    end
    return;
end

if string(scale) ~= "auto"
    error('plotPrtSystem3D:InvalidPolarizationScale', ...
        'PolarizationScale must be "auto" or a positive scalar.');
end

points = plottedExtentPoints(T, outputs, options);
if isempty(points)
    scale = 0.2;
    return;
end

span = max(max(points, [], 1) - min(points, [], 1));
scale = max(options.PolarizationScaleFraction * span, eps);
end

function points = plottedExtentPoints(T, outputs, options)
points = zeros(0,3);
if options.ShowSystem
    numSurfaces = height(T) - 1;
    vertexZ = prtSurfaceVertexZ(T);
    for surfaceIndex = 1:numSurfaces
        clearAperture = prtClearAperture(T(surfaceIndex,:), options.DefaultClearAperture);
        z = vertexZ(surfaceIndex);
        points = [points; ...
            -clearAperture, -clearAperture, z; ...
             clearAperture,  clearAperture, z];
    end
end

for outputIndex = 1:numel(outputs)
    rayOutput = outputs{outputIndex};
    for finalRayId = rayOutput.finalRayIds(:).'
        history = rayOutput.rays(finalRayId).history;
        points = [points; historyPositions(rayOutput, history)];
    end
end
end

function units = prtLengthUnits(T)
units = "";
if isfield(T.Properties.UserData, 'lambdaUnits')
    units = string(T.Properties.UserData.lambdaUnits);
end
if strlength(units) == 0
    units = "length units";
end
end

function padAxes3(points, paddingFraction)
mins = min(points, [], 1);
maxs = max(points, [], 1);
span = max(maxs - mins, eps);
padding = paddingFraction * span;
xlim([mins(1) - padding(1), maxs(1) + padding(1)]);
ylim([mins(2) - padding(2), maxs(2) + padding(2)]);
zlim([mins(3) - padding(3), maxs(3) + padding(3)]);
end

function addModeLegend(options)
legendHandles = gobjects(0);
legendLabels = strings(0);
modeNames = ["input", "ordinary", "extraordinary", "reflected"];
for modeName = modeNames
    if isfield(options.ModeColors, char(modeName))
        color = options.ModeColors.(char(modeName));
        legendHandles(end+1) = plot3(nan, nan, nan, '-', ...
            'Color', color, 'LineWidth', options.LineWidth); %#ok<AGROW>
        legendLabels(end+1) = modeName; %#ok<AGROW>
    end
end
legendHandles(end+1) = plot3(nan, nan, nan, '-', ...
    'Color', options.PolarizationColor, 'LineWidth', options.GlyphLineWidth);
legendLabels(end+1) = "polarization";
legend(legendHandles, legendLabels, 'Location', 'best');
end
