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

rayTraceData = prtInitializeRayTraceData(T, options);
rayTraceData.rays = prtMakeInitialRay(k_inc, x_inc, Ein, ...
    string(T.MaterialType(1)), T.IndexData{1}, struct('currentVertexZ', 0));

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
        ray = prtPropagateRayToHit( ...
            ray, hit, T.Properties.UserData.lambda, false, ...
            options.encodePropagationPhaseInP);
        rayTraceData.rays(rayId) = ray;

        interaction = traceSurfaceInteraction( ...
            T, surfaceIndex, ray, hit, normal, options);

        interaction.interceptData = interceptData;
        interaction.parentRayId = rayId;
        interaction.surfaceIndex = surfaceIndex;
        interaction.children = attachCurrentVertexZ(interaction.children, interceptData.targetZ);
        rayTraceData.interactions(end+1) = interaction;

        childIds = prtAppendChildRays(rayTraceData.rays, rayId, interaction, options);
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

function children = attachCurrentVertexZ(children, targetZ)
for childIndex = 1:numel(children)
    children(childIndex).metadata.currentVertexZ = targetZ;
end
end
