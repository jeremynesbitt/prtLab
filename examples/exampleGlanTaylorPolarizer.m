% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

%EXAMPLEGLANTAYLORPOLARIZER Trace a calcite Glan-Taylor polarizer.
%
% The air-gap incidence angle lies between the ordinary and extraordinary
% critical angles. The ordinary branch therefore undergoes TIR while the
% extraordinary branch is transmitted into the second prism.
exampleDir = fileparts(mfilename('fullpath'));
addpath(fullfile(exampleDir, '..', 'prtLab'));

lambda = 0.633; % um
calcite = struct('nO', 1.656, 'nE', 1.485);
opticAxis = [1;0;0];

gapIncidenceAngle = 40; % deg
legY = 10;              % um
legZ = legY*tand(gapIncidenceAngle);

scene = createGlanTaylorPolarizerScene(calcite, ...
    LegY=legY, ...
    LegZ=legZ, ...
    Width=8, ...
    AirGap=0.2, ...
    OpticAxis=opticAxis, ...
    lambda=lambda, ...
    lambdaUnits="um", ...
    Name="Calcite Glan-Taylor polarizer");

thetaCriticalO = asind(scene.outside.IndexData.n/calcite.nO);
thetaCriticalE = asind(scene.outside.IndexData.n/calcite.nE);
fprintf('Gap incidence angle: %.3f deg\n', gapIncidenceAngle);
fprintf('Ordinary critical angle: %.3f deg\n', thetaCriticalO);
fprintf('Extraordinary critical angle: %.3f deg\n', thetaCriticalE);

kIn = [0;0;1];
xIn = [0;legY/2;-1];
Ein = [1;1]/sqrt(2); % Excite both crystal eigenmodes.

options = prtDefaultOptions();
options.maxInteractions = 8;
% Keep the two dominant polarization paths while suppressing weak ghosts.
options.minAmplitude = 0.2;
options.minRelativeFlux = 1e-3;

rayOutput = polarizationRayTraceScene( ...
    scene, kIn, xIn, Ein, options);

disp("Hit sequence:");
for ii = 1:numel(rayOutput.interactions)
    interaction = rayOutput.interactions(ii);
    data = interaction.interceptData;
    parent = rayOutput.rays(interaction.parentRayId);
    fprintf('  %-24s / %-10s (%-8s), incident mode: %s\n', ...
        data.solidName, data.faceName, data.direction, parent.mode);
end

disp("Final branches:");
for rayId = rayOutput.finalRayIds
    ray = rayOutput.rays(rayId);
    fprintf('  ray %d: %-13s %-11s flux = %.6g\n', ...
        rayId, ray.mode, ray.branchType, ray.flux);
end

figure;
plotPrtSystem3D(scene, rayOutput, ...
    PreExtend=1, ...
    PostExtend=1, ...
    ShowPolarization=true, ...
    AxisPaddingFraction=0.2);
title(scene.name);
