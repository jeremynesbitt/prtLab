% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

%EXAMPLEFRESNELRHOMBSOLID Trace a ray through a simple Fresnel rhomb solid.

addpath(genpath('..'));

lambda = 0.589; % um
nGlass = 1.495;

rhomb = createFresnelRhombSolid(nGlass, ...
    Length=.5*28, ...
    Height=8, ...
    Width=8, ...
    Shear=.5*22.0417, ...
    lambda=lambda, ...
    lambdaUnits="um");
tta = .1;
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

options = prtDefaultOptions();
options.maxInteractions = 8;
options.minAmplitude = 1e-4;

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

disp("Hit sequence:");
for ii = 1:numel(rayOutput.interactions)
    data = rayOutput.interactions(ii).interceptData;
    fprintf("  %d: %s (%s)\n", ii, data.faceName, data.direction);
end

figure;
plotPrtSystem3D(rhomb, rayOutput, ...
    PostExtend=4, ...
    AxisPaddingFraction=0.22, ...
    SurfaceAlpha=0.28, ...
    ShowLegend=true);
title('Fresnel rhomb solid trace');
