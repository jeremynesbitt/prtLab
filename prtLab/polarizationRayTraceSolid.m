function rayTraceData = polarizationRayTraceSolid(solid, k_inc, x_inc, Ein, options)
%POLARIZATIONRAYTRACESOLID Trace a ray through an isotropic polygon solid.
%
%   This first solid tracer supports isotropic media and planar polygon
%   faces. Face ordering is geometric: the nearest forward face is used at
%   each step.

arguments
    solid (1,1) struct
    k_inc (3,1) double
    x_inc (3,1) double
    Ein double
    options struct = struct()
end

options = solidDefaultOptions(options);
k_inc = k_inc / norm(k_inc);

rayTraceData = prtInitializeRayTraceData(solid, options);
initialMetadata = struct('mediumName', "outside", 'lambda', solid.lambda);
rayTraceData.rays = prtMakeInitialRay(k_inc, x_inc, Ein, ...
    string(solid.outside.MaterialType), solid.outside.IndexData, initialMetadata);

activeIds = 1;
for bounceIndex = 1:options.maxInteractions
    nextActiveIds = [];
    for ii = 1:numel(activeIds)
        rayId = activeIds(ii);
        ray = rayTraceData.rays(rayId);
        if ~ray.active
            continue;
        end

        [hit, faceIndex, distance] = intersectSolidFaces(solid, ray, options);
        if isempty(faceIndex)
            rayTraceData.rays(rayId).active = true;
            nextActiveIds(end+1) = rayId; %#ok<AGROW>
            continue;
        end

        ray = prtPropagateRayToHit(ray, hit, solid.lambda, true, true);
        rayTraceData.rays(rayId) = ray;

        face = solid.faces(faceIndex);
        [localSystem, surfaceNormal, directionName] = localInterfaceSystem(solid, ray, face);
        interaction = traceSurfaceInteraction(localSystem, 1, ray, hit, surfaceNormal, options);
        interaction.parentRayId = rayId;
        interaction.surfaceIndex = faceIndex;
        interaction.interceptData = struct( ...
            'faceIndex', faceIndex, ...
            'faceName', face.name, ...
            'distance', distance, ...
            'direction', directionName);

        rayTraceData.interactions(end+1) = interaction;
        childMetadata = struct( ...
            'solidFaceIndex', interaction.surfaceIndex, ...
            'solidFaceName', interaction.interceptData.faceName);
        childResult = prtAppendChildRays( ...
            rayTraceData.rays, rayId, interaction, options, childMetadata);
        rayTraceData.rays = childResult.rays;
        nextActiveIds = [nextActiveIds, childResult.ids]; %#ok<AGROW>
        rayTraceData.rays(rayId).active = false;
    end

    activeIds = nextActiveIds;
    if isempty(activeIds)
        break;
    end
end

rayTraceData.finalRayIds = activeIds;

end

function options = solidDefaultOptions(userOptions)
options = prtDefaultOptions(userOptions);
if ~isfield(options, 'maxInteractions')
    options.maxInteractions = 12;
end
if ~isfield(options, 'faceTolerance')
    options.faceTolerance = 1e-9;
end
end

function [hit, faceIndex, distance] = intersectSolidFaces(solid, ray, options)
hit = [];
faceIndex = [];
distance = inf;

for ii = 1:numel(solid.faces)
    face = solid.faces(ii);
    denom = dot(ray.k, face.normal);
    if abs(denom) < options.faceTolerance
        continue;
    end
    t = dot(face.point - ray.position, face.normal) / denom;
    if t <= options.faceTolerance || t >= distance
        continue;
    end
    candidate = ray.position + t*ray.k;
    if isPointInFace(candidate, face, options.faceTolerance)
        hit = candidate;
        faceIndex = ii;
        distance = t;
    end
end
end

function tf = isPointInFace(point, face, tol)
vertices = face.vertices;
normal = face.normal;
tf = true;
numVertices = size(vertices, 2);
for jj = 1:numVertices
    v1 = vertices(:,jj);
    v2 = vertices(:,mod(jj, numVertices) + 1);
    if dot(cross(v2 - v1, point - v1), normal) < -tol
        tf = false;
        return;
    end
end
end

function [T, surfaceNormal, directionName] = localInterfaceSystem(solid, ray, face)
entering = dot(ray.k, face.normal) < 0;
if entering
    mediumIn = solid.outside;
    mediumOut = solid.material;
    surfaceNormal = -face.normal;
    directionName = "entering";
else
    mediumIn = solid.material;
    mediumOut = solid.outside;
    surfaceNormal = face.normal;
    directionName = "exiting";
end

T = createOpticalSystem(solid.lambda);
T.Properties.UserData.lambdaUnits = solid.lambdaUnits;
T = addSurface(T, Inf, 0, "plane", faceSurfaceData(face), ...
    mediumIn.MaterialType, mediumIn.IndexData, mediumIn.AxisData, face.CoatingData);
T = addSurface(T, Inf, 0, "plane", struct(), ...
    mediumOut.MaterialType, mediumOut.IndexData, mediumOut.AxisData, struct('type', 'bare'));
end

function surfaceData = faceSurfaceData(face)
surfaceData = struct();
if isfield(face, 'SurfaceData')
    surfaceData = face.SurfaceData;
end
end
