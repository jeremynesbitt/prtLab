function rayTraceData = polarizationRayTrace(T, k_inc, x_inc, Ein, options)
%POLARIZATIONRAYTRACE Trace a ray through a prtLab optical system.
%
%   rayTraceData = polarizationRayTrace(T, k_inc, x_inc, Ein)
%   rayTraceData = polarizationRayTrace(T, k_inc, x_inc, Ein, options)
%
%
%   Ein may be a 2x1 Jones vector or a 3x1 electric field. A Jones vector is
%   interpreted in the input double-pole basis stored in ray.localBasis.

arguments
    T table
    k_inc (3,1) double
    x_inc (3,1) double
    Ein double
    options struct = prtDefaultOptions()
end

options = prtDefaultOptions(options);
k_inc = k_inc / norm(k_inc);

rayTraceData = initializeRayTraceData(T, options);
rayTraceData.rays = makeInitialRay(k_inc, x_inc, Ein, T);

activeIds = 1;
for surfaceIndex = 1:height(T)-1
    nextActiveIds = [];

    for ii = 1:numel(activeIds)
        rayId = activeIds(ii);
        ray = rayTraceData.rays(rayId);

        if ~ray.active
            continue;
        end

        [hit, normal, interceptData] = intersectSurface(T(surfaceIndex,:), ray);

        interaction = traceSurfaceInteraction( ...
            T, surfaceIndex, ray, hit, normal, options);

        interaction.interceptData = interceptData;
        interaction.parentRayId = rayId;
        interaction.surfaceIndex = surfaceIndex;
        interaction.children = attachCurrentVertexZ(interaction.children, interceptData.targetZ);
        rayTraceData.interactions(end+1) = interaction;

        childIds = appendChildRays(rayTraceData.rays, rayId, interaction, options);
        rayTraceData.rays = childIds.rays;
        nextActiveIds = [nextActiveIds, childIds.ids]; %#ok<AGROW>

        rayTraceData.rays(rayId).active = false;
    end

    activeIds = nextActiveIds;
    if isempty(activeIds)
        break;
    end
end

rayTraceData.finalRayIds = activeIds;

end

function rayTraceData = initializeRayTraceData(T, options)
rayTraceData = struct();
rayTraceData.system = T;
rayTraceData.options = options;
rayTraceData.rays = repmat(emptyRayBranch(), 0, 1);
rayTraceData.interactions = repmat(emptyInteractionRecord(), 0, 1);
rayTraceData.finalRayIds = [];
end

function ray = makeInitialRay(k_inc, x_inc, Ein, T)
ray = emptyRayBranch();
ray.id = 1;
ray.parentId = 0;
ray.surfaceIndex = 0;
ray.mode = "input";
ray.mediumType = string(T.MaterialType(1));
ray.position = x_inc;
ray.k = k_inc;
ray.S = k_inc;
[fieldE, localBasis] = initializeInputElectricField(Ein, k_inc);
nInc = incidentIndex(T(1,:));
modeE = fieldE / norm(fieldE);
modeH = nInc * makeK(k_inc) * modeE;
fieldH = nInc * makeK(k_inc) * fieldE;
ray.modeE = modeE;
ray.modeH = modeH;
ray.fieldE = fieldE;
ray.fieldH = fieldH;
ray.E = fieldE;
ray.H = fieldH;
ray.P = eye(3);
ray.Q = eye(3);
ray.O = eye(3);
ray.localBasis = localBasis;
ray.amplitude = norm(fieldE);
ray.flux = real(dot(cross(fieldE, conj(fieldH)), k_inc));
ray.OPL = 0;
ray.active = true;
ray.history = 1;
ray.metadata = struct('n', nInc, 'currentVertexZ', 0);
end

function [fieldE, localBasis] = initializeInputElectricField(Ein, k_inc)
if numel(Ein) == 2
    Ein = Ein(:);
    xPole = [1;0;0];
    [x_dp, y_dp] = doublePoleBasisVectors(k_inc, [0;0;1], xPole);
    x_dp = x_dp(:);
    y_dp = y_dp(:);
    fieldE = Ein(1)*x_dp + Ein(2)*y_dp;
    localBasis = struct( ...
        'x', x_dp, ...
        'y', y_dp, ...
        'basisDirection', k_inc, ...
        'inputConvention', "doublePoleJones");
elseif numel(Ein) == 3
    fieldE = Ein(:);
    longitudinal = dot(fieldE, k_inc);
    if abs(longitudinal) > 100*eps(max(1,norm(fieldE)))
        fieldE = fieldE - longitudinal*k_inc;
    end
    fieldE = fieldE(:);
    localBasis = struct( ...
        'x', [], ...
        'y', [], ...
        'basisDirection', k_inc, ...
        'inputConvention', "global3D");
else
    error('polarizationRayTrace:InvalidInputField', ...
        'Ein must be either a 2x1 Jones vector or a 3x1 electric field.');
end

if norm(fieldE) == 0
    error('polarizationRayTrace:ZeroInputField', ...
        'Ein must have nonzero electric-field amplitude.');
end
end

function nInc = incidentIndex(surfaceRow)
indexData = surfaceRow.IndexData{1};
if isfield(indexData, 'n')
    nInc = indexData.n;
elseif isfield(indexData, 'nO')
    nInc = indexData.nO;
else
    error('polarizationRayTrace:MissingIndex', ...
        'The first surface IndexData must contain n or nO.');
end
end

function childResult = appendChildRays(rays, parentRayId, interaction, options)
childIds = [];

for ii = 1:numel(interaction.children)
    child = interaction.children(ii);

    if child.flux < options.minFlux
        continue;
    end
    if abs(child.amplitude) < options.minAmplitude
        continue;
    end

    child.id = numel(rays) + 1;
    child.parentId = parentRayId;
    child.history = [rays(parentRayId).history, child.id];
    rays(end+1) = child; %#ok<AGROW>
    childIds(end+1) = child.id; %#ok<AGROW>
end

childResult = struct('rays', rays, 'ids', childIds);
end

function children = attachCurrentVertexZ(children, targetZ)
for childIndex = 1:numel(children)
    children(childIndex).metadata.currentVertexZ = targetZ;
end
end
