% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function tests = rayBranchPolicyTest
%RAYBRANCHPOLICYTEST Branch generation and tracer pruning tests.

tests = functiontests(localfunctions);

end

function testSequentialTraceGeneratesPhysicalReflections(testCase)
T = simpleGlassPlate();
output = polarizationRayTrace(T, [0;0;1], [0;0;-1], [1;0]);

verifyGreaterThanOrEqual(testCase, numel(output.interactions), 2);
verifyEqual(testCase, [output.interactions(1).children.branchType], ...
    ["transmitted", "reflected"]);
verifyTrue(testCase, any([output.rays.branchType] == "reflected"));
end

function testRelativeFluxPrunesWeakFresnelReflections(testCase)
T = simpleGlassPlate();
options = prtDefaultOptions(struct('minRelativeFlux', 0.1));
output = polarizationRayTrace( ...
    T, [0;0;1], [0;0;-1], [1;0], options);

verifyEqual(testCase, numel(output.finalRayIds), 1);
verifyEqual(testCase, ...
    output.rays(output.finalRayIds).branchType, "transmitted");
verifyEqual(testCase, [output.interactions(1).children.branchType], ...
    ["transmitted", "reflected"]);
end

function testMaxBranchesCapsRegisteredRayTree(testCase)
T = simpleGlassPlate();
options = prtDefaultOptions(struct('maxBranches', 2));
output = polarizationRayTrace( ...
    T, [0;0;1], [0;0;-1], [1;0], options);

verifyLessThanOrEqual(testCase, numel(output.rays), 2);
end

function T = simpleGlassPlate()
T = createOpticalSystem(0.55);
air = struct('n', 1.0);
glass = struct('n', 1.5);
coating = struct('type', 'bare');

T = addSurface(T, Inf, 0, "plane", struct(), ...
    "isotropic", air, struct(), coating);
T = addSurface(T, Inf, 1, "plane", struct(), ...
    "isotropic", glass, struct(), coating);
T = addSurface(T, Inf, 0, "plane", struct(), ...
    "isotropic", air, struct(), coating);
end
