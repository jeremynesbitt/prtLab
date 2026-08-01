% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function outputPath = generateSceneTraceBaselines(outputPath)
%GENERATESCENETRACEBASELINES Export branched solid/scene trace references.

repositoryRoot = fileparts(fileparts(mfilename('fullpath')));
if nargin < 1
    outputPath = fullfile(repositoryRoot, 'pyprtLab', 'tests', ...
        'baselines', 'sceneTraceBaselines.json');
end
addpath(fullfile(repositoryRoot, 'prtLab'));

fixture = struct();
fixture.schemaVersion = 1;
fixture.referenceImplementation = "MATLAB prtLab";
fixture.referenceCommit = "2682ac1deb00040ee59f48b8516a962a6b6a24ff";
fixture.referenceWorkingTreeChanges = [ ...
    "prtLab/uniaxialModesFromTangentialQ.m", ...
    "prtLab/traceUniaxialToUniaxial.m"];
fixture.generator = "tests/generateSceneTraceBaselines.m";
fixture.absoluteTolerance = 5e-10;
fixture.fresnelRhomb = fresnelRhombCases();
fixture.glanTaylor = glanTaylorCase();

outputFolder = fileparts(outputPath);
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end
fileId = fopen(outputPath, 'w');
cleanup = onCleanup(@() fclose(fileId));
fprintf(fileId, '%s\n', jsonencode(fixture, PrettyPrint=true));
fprintf('Wrote %s\n', outputPath);
end


function cases = fresnelRhombCases()
lambda = 0.589;
nGlass = 1.4965;
rhomb = createFresnelRhombSolid(nGlass, ...
    Length=14, Height=8, Width=8, Shear=11.02085, ...
    lambda=lambda, lambdaUnits="um");
options = prtDefaultOptions();
options.maxInteractions = 8;
options.minAmplitude = 1e-4;
coord = struct('type', "doublePole", ...
    'a_loc', [0;0;1], 'x_o', [1;0;0]);

angles = [0, 5];
for index = numel(angles):-1:1
    angle = angles(index);
    kIn = [0; sind(angle); cosd(angle)];
    output = polarizationRayTraceSolid( ...
        rhomb, kIn, [0;5;-2], [1;1]/sqrt(2), options);
    [ray, rayId] = selectDominantFinalRay(output);
    jones = transformPtoJones(ray.P, kIn, ray.k, coord);
    cases(index) = struct( ...
        'angleDeg', angle, ...
        'rayCount', numel(output.rays), ...
        'interactionCount', numel(output.interactions), ...
        'finalRayCount', numel(output.finalRayIds), ...
        'dominantRayId', rayId, ...
        'dominant', serializeRay(output, ray), ...
        'jones', packNumeric(jones), ...
        'retardance', packNumeric(calcRetardanceFromJonesMatrix(jones)));
end
end


function result = glanTaylorCase()
lambda = 0.633;
indexData = struct('nO', 1.656, 'nE', 1.485);
legY = 10;
scene = createGlanTaylorPolarizerScene(indexData, ...
    LegY=legY, LegZ=legY*tand(40), Width=8, AirGap=0.2, ...
    OpticAxis=[1;0;0], lambda=lambda, lambdaUnits="um", ...
    Name="Calcite Glan-Taylor polarizer");
options = prtDefaultOptions();
options.maxInteractions = 8;
options.minAmplitude = 0.2;
options.minRelativeFlux = 1e-3;
output = polarizationRayTraceScene( ...
    scene, [0;0;1], [0;legY/2;-1], [1;1]/sqrt(2), options);

finalRays = serializeRay(output, output.rays(output.finalRayIds(1)));
for index = 2:numel(output.finalRayIds)
    finalRays(end+1) = serializeRay( ...
        output, output.rays(output.finalRayIds(index))); %#ok<AGROW>
end
result = struct( ...
    'rayCount', numel(output.rays), ...
    'interactionCount', numel(output.interactions), ...
    'finalRayIds', output.finalRayIds, ...
    'finalRays', finalRays);
end


function record = serializeRay(output, ray)
[faces, directions, cases, incidentModes] = rayPath(output, ray);
record = struct( ...
    'id', ray.id, ...
    'mode', ray.mode, ...
    'branchType', ray.branchType, ...
    'history', ray.history, ...
    'faces', faces, ...
    'directions', directions, ...
    'interfaceCases', cases, ...
    'incidentModes', incidentModes, ...
    'position', packNumeric(ray.position), ...
    'k', packNumeric(ray.k), ...
    'S', packNumeric(ray.S), ...
    'fieldE', packNumeric(ray.fieldE), ...
    'fieldH', packNumeric(ray.fieldH), ...
    'P', packNumeric(ray.P), ...
    'Q', packNumeric(ray.Q), ...
    'amplitude', packNumeric(ray.amplitude), ...
    'flux', packNumeric(ray.flux), ...
    'OPL', packNumeric(ray.OPL));
end


function [faces, directions, cases, modes] = rayPath(output, ray)
faces = strings(1, 0);
directions = strings(1, 0);
cases = strings(1, 0);
modes = strings(1, 0);
for rayId = ray.history
    index = find([output.interactions.parentRayId] == rayId, 1);
    if isempty(index)
        continue;
    end
    interaction = output.interactions(index);
    faces(end+1) = interaction.interceptData.faceName; %#ok<AGROW>
    directions(end+1) = interaction.interceptData.direction; %#ok<AGROW>
    cases(end+1) = interaction.caseName; %#ok<AGROW>
    modes(end+1) = interaction.incident.mode; %#ok<AGROW>
end
end


function packed = packNumeric(value)
packed = struct('real', real(value), 'imag', imag(value));
end
