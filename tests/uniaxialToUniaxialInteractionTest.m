% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function tests = uniaxialToUniaxialInteractionTest
%UNIAXIALTOUNIAXIALINTERACTIONTEST Verify oblique anisotropic boundaries.

tests = functiontests(localfunctions);
end


function testObliqueTangentialWavevectorAndFlux(testCase)
[T, entry, options] = makeObliqueCase();
normal = [0;0;1];

for incidentIndex = 1:2
    ray = entry.children(incidentIndex);
    interaction = traceSurfaceInteraction( ...
        T, 2, ray, [0;0;1], normal, options);

    verifyEqual(testCase, numel(interaction.children), 4);
    verifyLessThan(testCase, ...
        norm(interaction.diagnostics.boundaryResidual_m), 1e-12);
    verifyEqual(testCase, ...
        real(interaction.diagnostics.totalFluxRatioField), 1, ...
        AbsTol=5e-12);

    qIncident = ray.metadata.n*ray.k;
    qTangential = qIncident - dot(qIncident, normal)*normal;
    for childIndex = 1:numel(interaction.children)
        child = interaction.children(childIndex);
        qChild = child.metadata.n*child.k;
        childTangential = qChild - dot(qChild, normal)*normal;
        verifyEqual(testCase, childTangential, qTangential, ...
            AbsTol=5e-12);
    end

    ordinary = interaction.children(1);
    verifyEqual(testCase, ordinary.S, ordinary.k, AbsTol=5e-12);
end
end


function [T, entry, options] = makeObliqueCase()
air = struct('n', 1.0);
index1 = struct('nO', 1.656, 'nE', 1.485);
index2 = struct('nO', 1.603, 'nE', 1.497);
axis1 = prtNorm([0.35;-0.60;0.72]);
axis2 = prtNorm([-0.50;0.40;0.77]);
bare = struct('type', 'bare');

T = createOpticalSystem(0.633);
T = addSurface(T, Inf, 0, "plane", struct(), "isotropic", ...
    air, struct(), bare);
T = addSurface(T, Inf, 1, "plane", struct(), "uniaxial", ...
    index1, struct('opticAxis', axis1), bare);
T = addSurface(T, Inf, 0, "plane", struct(), "uniaxial", ...
    index2, struct('opticAxis', axis2), bare);

options = prtDefaultOptions();
initial = prtMakeInitialRay(calcKFromThetaXThetaY(20, 10), ...
    [0;0;0], [1;1]/sqrt(2), "isotropic", air);
entry = traceSurfaceInteraction(T, 1, initial, ...
    [0;0;0], [0;0;1], options);
end
