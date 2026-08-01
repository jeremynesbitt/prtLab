% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function rayTraceData = polarizationRayTraceSolid(solid, k_inc, x_inc, Ein, options)
%POLARIZATIONRAYTRACESOLID Trace a ray through one homogeneous solid.
%
%   This compatibility wrapper delegates to polarizationRayTraceScene.

arguments
    solid (1,1) struct
    k_inc (3,1) double
    x_inc (3,1) double
    Ein double
    options struct = struct()
end

if ~isfield(options, 'maxInteractions')
    options.maxInteractions = 12;
end
scene = createPrtScene(solid.lambda, solid.outside, ...
    lambdaUnits=solid.lambdaUnits, Name=solid.name);
scene = addSolidToPrtScene(scene, solid);
rayTraceData = polarizationRayTraceScene(scene, k_inc, x_inc, Ein, options);

% Preserve the historical single-solid system field for callers.
rayTraceData.system = solid;

end
