function tests = multiSolidSceneTest
%MULTISOLIDSCENETEST Physics and compatibility tests for scene tracing.

tests = functiontests(localfunctions);

end

function testTwoPrismAirGap(testCase)
gap = 0.2;
scene = createGlanTaylorPolarizerScene(struct('n', 1.2), ...
    LegY=10, LegZ=10, Width=8, AirGap=gap, lambda=0.633);

h1 = namedFace(scene.solids{1}, "hypotenuse");
h2 = namedFace(scene.solids{2}, "hypotenuse");
verifyEqual(testCase, ...
    abs(dot(h2.point - h1.point, h1.normal)), gap, 'AbsTol', 1e-12);

options = prtDefaultOptions();
options.maxInteractions = 8;
options.minAmplitude = 1e-8;
options.minRelativeFlux = 0.02;
output = polarizationRayTraceScene( ...
    scene, [0;0;1], [0;5;-1], [1;0], options);

finalRay = allTransmittedFinalRay(output);
verifyEqual(testCase, finalRay.metadata.solidIndex, 0);
verifyEqual(testCase, finalRay.metadata.mediumName, "outside");
verifyTrue(testCase, isfinite(finalRay.OPL));

history = finalRay.history;
parentIds = history(1:end-1);
interactionIndices = arrayfun(@(id) ...
    find([output.interactions.parentRayId] == id, 1), parentIds);
data = [output.interactions(interactionIndices).interceptData];
verifyEqual(testCase, string({data.faceName}), ...
    ["z leg", "hypotenuse", "hypotenuse", "z leg"]);
verifyEqual(testCase, [data.solidIndex], [1 1 2 2]);
verifyEqual(testCase, string({data.direction}), ...
    ["entering", "exiting", "entering", "exiting"]);

verifyEqual(testCase, [output.rays(history).mode], ...
    repmat("isotropic", 1, numel(history)));
verifyEqual(testCase, output.rays(history(1)).branchType, "input");
verifyEqual(testCase, [output.rays(history(2:end)).branchType], ...
    repmat("transmitted", 1, numel(history)-1));
end

function ray = allTransmittedFinalRay(output)
for rayId = output.finalRayIds
    candidate = output.rays(rayId);
    branchTypes = [output.rays(candidate.history).branchType];
    if all(branchTypes(2:end) == "transmitted")
        ray = candidate;
        return;
    end
end
error('multiSolidSceneTest:MissingTransmittedPath', ...
    'No all-transmitted final ray was found.');
end

function testSingleSolidWrapperMatchesScene(testCase)
solid = createFresnelRhombSolid(1.495, ...
    Length=14, Height=8, Width=8, Shear=11.02085, lambda=0.589);
options = prtDefaultOptions();
options.maxInteractions = 8;
options.minAmplitude = 1e-4;
kIn = [0;sin(0.1*pi/180);cos(0.1*pi/180)];
xIn = [0;5;-2];
Ein = [1;1]/sqrt(2);

wrapped = polarizationRayTraceSolid(solid, kIn, xIn, Ein, options);
scene = createPrtScene(solid.lambda, solid.outside, ...
    lambdaUnits=solid.lambdaUnits);
scene = addSolidToPrtScene(scene, solid);
direct = polarizationRayTraceScene(scene, kIn, xIn, Ein, options);

verifyEqual(testCase, numel(wrapped.interactions), numel(direct.interactions));
verifyEqual(testCase, wrapped.finalRayIds, direct.finalRayIds);
for ii = 1:numel(wrapped.rays)
    verifyEqual(testCase, wrapped.rays(ii).position, direct.rays(ii).position, ...
        'AbsTol', 1e-13);
    verifyEqual(testCase, wrapped.rays(ii).k, direct.rays(ii).k, ...
        'AbsTol', 1e-13);
    verifyEqual(testCase, wrapped.rays(ii).P, direct.rays(ii).P, ...
        'AbsTol', 1e-13);
    verifyEqual(testCase, wrapped.rays(ii).OPL, direct.rays(ii).OPL, ...
        'AbsTol', 1e-13);
end

allChildren = vertcat(wrapped.interactions.children);
reflected = allChildren([allChildren.branchType] == "reflected");
verifyNotEmpty(testCase, reflected);
verifyEqual(testCase, [reflected.mode], ...
    repmat("isotropic", 1, numel(reflected)));
end

function testRigidTransform(testCase)
solid = createRightTriangularPrismSolid(1.5, ...
    LegY=3, LegZ=4, Width=2, OpticAxis=[0;1;0]);
rotation = [1 0 0; 0 0 -1; 0 1 0];
translation = [2;3;4];
originalPoint = solid.faces(1).point;
originalNormal = solid.faces(1).normal;

transformed = transformPrtSolid(solid, rotation, translation);
verifyEqual(testCase, transformed.faces(1).point, ...
    rotation*originalPoint + translation, 'AbsTol', 1e-13);
verifyEqual(testCase, transformed.faces(1).normal, ...
    rotation*originalNormal, 'AbsTol', 1e-13);
end

function face = namedFace(solid, name)
index = find([solid.faces.name] == name, 1);
face = solid.faces(index);
end
