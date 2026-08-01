% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function tests = jonesAnalysisTest
%JONESANALYSISTEST Test Jones ellipse parameter calculations.

tests = functiontests(localfunctions);
end


function testLinearStates(testCase)
fields = [1,0,1/sqrt(2); 0,1,1/sqrt(2)];
[psi, ellipticity, ~, major, minor] = ...
    jonesMajorAxisEllipticity(fields);
verifyEqual(testCase, psi, [0,pi/2,pi/4], AbsTol=1e-14);
verifyEqual(testCase, ellipticity, zeros(1,3), AbsTol=1e-14);
verifyEqual(testCase, major, ones(1,3), AbsTol=1e-14);
verifyEqual(testCase, minor, zeros(1,3), AbsTol=1e-14);
end


function testCircularState(testCase)
[~, ellipticity, ~, major, minor] = ...
    jonesMajorAxisEllipticity([1;1i]/sqrt(2));
verifyEqual(testCase, ellipticity, 1, AbsTol=1e-14);
verifyEqual(testCase, major, minor, AbsTol=1e-14);
end
