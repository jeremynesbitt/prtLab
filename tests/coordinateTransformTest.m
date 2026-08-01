% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function tests = coordinateTransformTest
%COORDINATETRANSFORMTEST Test P-to-Jones coordinate conventions.

tests = functiontests(localfunctions);
end


function testDipolePoleReturnsOrthonormalBasis(testCase)
[xLocal, yLocal] = dipoleBasisVectors([0,0,1], [0,0,1]);
verifyEqual(testCase, norm(xLocal), 1, AbsTol=1e-14);
verifyEqual(testCase, norm(yLocal), 1, AbsTol=1e-14);
verifyEqual(testCase, dot(xLocal, yLocal), 0, AbsTol=1e-14);
verifyEqual(testCase, dot(xLocal, [0,0,1]), 0, AbsTol=1e-14);
verifyEqual(testCase, dot(yLocal, [0,0,1]), 0, AbsTol=1e-14);
end


function testIdentityPIsIdentityJonesForSharedCoordinates(testCase)
k = prtNorm([0.2;-0.3;0.93]);
coordinates = { ...
    struct('type', "doublePole", 'a_loc', [0;0;1], 'x_o', [1;0;0]), ...
    struct('type', "dipole", 'a_loc', [0;0;1]), ...
    struct('type', "sp", 'normal', [0;0;1])};
for index = 1:numel(coordinates)
    J = transformPtoJones(eye(3), k, k, coordinates{index});
    verifyEqual(testCase, J, eye(2), AbsTol=1e-14);
end
end


function testSPNormalIncidenceIsRejected(testCase)
coord = struct('type', "sp", 'normal', [0;0;1]);
verifyError(testCase, ...
    @() transformPtoJones(eye(3), [0;0;1], [0;0;1], coord), ...
    'transformPtoJones:SingularSPBasis');
end
