% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function scene = addSolidToPrtScene(scene, solid)
%ADDSOLIDTOPRTSCENE Add a homogeneous solid to a prtLab scene.

arguments
    scene (1,1) struct
    solid (1,1) struct
end

if ~isfield(scene, 'type') || string(scene.type) ~= "scene"
    error('addSolidToPrtScene:InvalidScene', ...
        'The first input must be a scene created by createPrtScene.');
end
if ~isfield(solid, 'type') || string(solid.type) ~= "solid" || ...
        ~all(isfield(solid, {'faces', 'material'}))
    error('addSolidToPrtScene:InvalidSolid', ...
        'The second input must be a prtLab solid.');
end
if abs(solid.lambda - scene.lambda) > 100*eps(max(scene.lambda, solid.lambda))
    error('addSolidToPrtScene:WavelengthMismatch', ...
        'All scene solids must use the scene wavelength.');
end
if string(solid.lambdaUnits) ~= string(scene.lambdaUnits)
    error('addSolidToPrtScene:UnitMismatch', ...
        'All scene solids must use the scene wavelength/length units.');
end

% The scene owns the shared medium between disjoint solids.
solid.outside = scene.outside;
scene.solids{end+1,1} = solid;

end
