% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function tests = propagationPhasePolicyTest
%PROPAGATIONPHASEPOLICYTEST OPL tracking and optional P-phase tests.

tests = functiontests(localfunctions);

end

function testPropagationPhaseIsOptIn(testCase)
T = exampleQuarterWavePlateSystem();
kInc = calcKFromThetaXThetaY(30, 30);
xInc = [0; 0; 0];
Ein = [0; 1];

defaultOutput = polarizationRayTrace(T, kInc, xInc, Ein);
phaseOutput = polarizationRayTrace(T, kInc, xInc, Ein, ...
    struct('encodePropagationPhaseInP', true));

verifyFalse(testCase, defaultOutput.options.encodePropagationPhaseInP);
verifyTrue(testCase, phaseOutput.options.encodePropagationPhaseInP);
verifyEqual(testCase, sortedFinalOPL(defaultOutput), ...
    sortedFinalOPL(phaseOutput), 'AbsTol', 1e-12);
verifyGreaterThan(testCase, ...
    norm(coherentFinalP(phaseOutput) - coherentFinalP(defaultOutput), 'fro'), ...
    1e-6);
end

function opl = sortedFinalOPL(output)
opl = sort([output.rays(output.finalRayIds).OPL]);
end

function P = coherentFinalP(output)
P = zeros(3);
for rayId = output.finalRayIds
    P = P + output.rays(rayId).P;
end
end
