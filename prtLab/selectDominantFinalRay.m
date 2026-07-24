function [ray, rayId] = selectDominantFinalRay(rayOutput)
%SELECTDOMINANTFINALRAY Return the final ray with the largest flux magnitude.

finalRayIds = rayOutput.finalRayIds;
if isempty(finalRayIds)
    error('selectDominantFinalRay:NoFinalRays', ...
        'The ray trace does not contain any final rays.');
end

flux = abs([rayOutput.rays(finalRayIds).flux]);
[~, index] = max(flux);
rayId = finalRayIds(index);
ray = rayOutput.rays(rayId);

end
