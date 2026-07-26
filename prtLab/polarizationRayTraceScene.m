function rayTraceData = polarizationRayTraceScene(scene, k_inc, x_inc, Ein, options)
%POLARIZATIONRAYTRACESCENE Trace a branched ray through multiple solids.
%
%   Input rays currently begin in scene.outside. Solids must be disjoint;
%   nested and overlapping solids are not supported.

arguments
    scene (1,1) struct
    k_inc (3,1) double
    x_inc (3,1) double
    Ein double
    options struct = struct()
end

validateScene(scene);
options = sceneDefaultOptions(options);
k_inc = k_inc / norm(k_inc);

rayTraceData = prtInitializeRayTraceData(scene, options);
initialMetadata = struct( ...
    'mediumName', "outside", ...
    'solidIndex', 0, ...
    'lambda', scene.lambda);
rayTraceData.rays = prtMakeInitialRay(k_inc, x_inc, Ein, ...
    string(scene.outside.MaterialType), scene.outside.IndexData, initialMetadata);

activeIds = 1;
escapedIds = [];
for interactionIndex = 1:options.maxInteractions
    nextActiveIds = [];
    for ii = 1:numel(activeIds)
        rayId = activeIds(ii);
        ray = rayTraceData.rays(rayId);
        if ~ray.active
            continue;
        end

        [hit, solidIndex, faceIndex, distance] = ...
            intersectScene(scene, ray, options);
        if isempty(solidIndex)
            escapedIds(end+1) = rayId; %#ok<AGROW>
            continue;
        end

        ray = prtPropagateRayToHit( ...
            ray, hit, scene.lambda, true, ...
            options.encodePropagationPhaseInP);
        rayTraceData.rays(rayId) = ray;

        solid = scene.solids{solidIndex};
        face = solid.faces(faceIndex);
        [localSystem, surfaceNormal, directionName, outputSolidIndex] = ...
            localInterfaceSystem(scene, solid, solidIndex, ray, face);
        interaction = traceSurfaceInteraction( ...
            localSystem, 1, ray, hit, surfaceNormal, options);
        interaction.parentRayId = rayId;
        interaction.surfaceIndex = faceIndex;
        interaction.interceptData = struct( ...
            'solidIndex', solidIndex, ...
            'solidName', solid.name, ...
            'faceIndex', faceIndex, ...
            'faceName', face.name, ...
            'distance', distance, ...
            'direction', directionName);
        rayTraceData.interactions(end+1) = interaction;

        childMetadata = struct( ...
            'solidFaceIndex', faceIndex, ...
            'solidFaceName', face.name, ...
            'sceneSolidIndex', solidIndex, ...
            'sceneSolidName', solid.name);
        childResult = prtAppendChildRays( ...
            rayTraceData.rays, rayId, interaction, options, childMetadata);
        rayTraceData.rays = childResult.rays;
        rayTraceData.rays = updateChildOwnership( ...
            rayTraceData.rays, childResult.ids, ray, ...
            outputSolidIndex, solid);
        nextActiveIds = [nextActiveIds, childResult.ids]; %#ok<AGROW>
        rayTraceData.rays(rayId).active = false;
    end

    activeIds = nextActiveIds;
    if isempty(activeIds)
        break;
    end
end

rayTraceData.finalRayIds = unique([escapedIds, activeIds], 'stable');

end

function options = sceneDefaultOptions(userOptions)
options = prtDefaultOptions(userOptions);
if ~isfield(options, 'maxInteractions')
    options.maxInteractions = 24;
end
if ~isfield(options, 'faceTolerance')
    options.faceTolerance = 1e-9;
end
end

function validateScene(scene)
if ~isfield(scene, 'type') || string(scene.type) ~= "scene" || ...
        ~all(isfield(scene, {'lambda', 'lambdaUnits', 'outside', 'solids'}))
    error('polarizationRayTraceScene:InvalidScene', ...
        'scene must be created by createPrtScene.');
end
end

function [hit, solidIndex, faceIndex, distance] = intersectScene(scene, ray, options)
hit = [];
solidIndex = [];
faceIndex = [];
distance = inf;
currentSolidIndex = raySolidIndex(ray);

if currentSolidIndex == 0
    solidIndices = 1:numel(scene.solids);
else
    solidIndices = currentSolidIndex;
end

for candidateSolidIndex = solidIndices
    solid = scene.solids{candidateSolidIndex};
    for candidateFaceIndex = 1:numel(solid.faces)
        face = solid.faces(candidateFaceIndex);
        denom = dot(ray.k, face.normal);
        if currentSolidIndex == 0
            if denom >= -options.faceTolerance
                continue;
            end
        elseif denom <= options.faceTolerance
            continue;
        end

        t = dot(face.point - ray.position, face.normal) / denom;
        if ~isreal(t) || t <= options.faceTolerance || t >= distance
            continue;
        end
        candidate = ray.position + t*ray.k;
        if isPointInFace(candidate, face, options.faceTolerance)
            hit = candidate;
            solidIndex = candidateSolidIndex;
            faceIndex = candidateFaceIndex;
            distance = t;
        end
    end
end
end

function tf = isPointInFace(point, face, tolerance)
vertices = face.vertices;
normal = face.normal;
tf = true;
numVertices = size(vertices, 2);
for ii = 1:numVertices
    v1 = vertices(:,ii);
    v2 = vertices(:,mod(ii, numVertices) + 1);
    if dot(cross(v2 - v1, point - v1), normal) < -tolerance
        tf = false;
        return;
    end
end
end

function [T, normal, directionName, outputSolidIndex] = ...
        localInterfaceSystem(scene, solid, solidIndex, ray, face)
currentSolidIndex = raySolidIndex(ray);
if currentSolidIndex == 0
    mediumIn = scene.outside;
    mediumOut = solid.material;
    normal = -face.normal;
    directionName = "entering";
    outputSolidIndex = solidIndex;
else
    mediumIn = solid.material;
    mediumOut = scene.outside;
    normal = face.normal;
    directionName = "exiting";
    outputSolidIndex = 0;
end

T = createOpticalSystem(scene.lambda);
T.Properties.UserData.lambdaUnits = scene.lambdaUnits;
T = addSurface(T, Inf, 0, "plane", faceSurfaceData(face), ...
    mediumIn.MaterialType, mediumIn.IndexData, mediumIn.AxisData, face.CoatingData);
T = addSurface(T, Inf, 0, "plane", struct(), ...
    mediumOut.MaterialType, mediumOut.IndexData, mediumOut.AxisData, ...
    struct('type', 'bare'));
end

function surfaceData = faceSurfaceData(face)
surfaceData = struct();
if isfield(face, 'SurfaceData')
    surfaceData = face.SurfaceData;
end
end

function rays = updateChildOwnership(rays, childIds, parentRay, ...
        transmittedSolidIndex, solid)
parentSolidIndex = raySolidIndex(parentRay);
for childId = childIds
    if rays(childId).branchType == "reflected"
        childSolidIndex = parentSolidIndex;
    else
        childSolidIndex = transmittedSolidIndex;
    end

    rays(childId).metadata.solidIndex = childSolidIndex;
    if childSolidIndex == 0
        rays(childId).metadata.mediumName = "outside";
    else
        rays(childId).metadata.mediumName = solid.name;
    end
end
end

function solidIndex = raySolidIndex(ray)
solidIndex = 0;
if isfield(ray.metadata, 'solidIndex') && ~isempty(ray.metadata.solidIndex)
    solidIndex = ray.metadata.solidIndex;
end
end
