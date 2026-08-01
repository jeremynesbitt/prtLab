% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function solid = createFresnelRhombSolid(n, options)
%CREATEFRESNELRHOMBSOLID Create a simple isotropic Fresnel rhomb solid.
%
%   solid = createFresnelRhombSolid(n) returns a parallelogram rhomb
%   extruded in x. Length units follow options.lambda.

arguments
    n (1,1) double {mustBePositive}
    options.Length (1,1) double {mustBePositive} = 20
    options.Height (1,1) double {mustBePositive} = 8
    options.Width (1,1) double {mustBePositive} = 8
    options.Shear (1,1) double = 5
    options.lambda (1,1) double {mustBePositive} = 0.633
    options.lambdaUnits (1,1) string = "um"
end

L = options.Length;
H = options.Height;
W = options.Width;
s = options.Shear;

% Cross-section in y-z. The first coordinate is y, second is z.
A = [0; 0];
B = [H; 0];
C = [H + s; L];
D = [s; L];
solid = prtCreateExtrudedPolygonSolid( ...
    "Fresnel rhomb", [A, B, C, D], ...
    ["entrance"; "top"; "exit"; "bottom"], W, ...
    isotropicMedium(1.0), isotropicMedium(n), ...
    options.lambda, options.lambdaUnits);

end

function medium = isotropicMedium(n)
medium = struct( ...
    'MaterialType', "isotropic", ...
    'IndexData', struct('n', n), ...
    'AxisData', struct());
end
