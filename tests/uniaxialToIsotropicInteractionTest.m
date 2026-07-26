function tests = uniaxialToIsotropicInteractionTest
%UNIAXIALTOISOTROPICINTERACTIONTEST Test transmitted and reflected branches.

tests = functiontests(localfunctions);

end


function testThreeBranchInteractionRecord(testCase)
T = exampleQuarterWavePlateSystem();
output = polarizationRayTrace(T, [0;0;1], [0;0;0], [0;1]);
interaction = firstInteraction(output, "uniaxialToIsotropic");

verifyEqual(testCase, numel(interaction.children), 3);
verifyEqual(testCase, [interaction.children.mode], ...
    ["isotropic", "ordinary", "extraordinary"]);
verifyEqual(testCase, [interaction.children.branchType], ...
    ["transmitted", "reflected", "reflected"]);
verifyTrue(testCase, all(isfield(interaction.P, ...
    {'transmitted', 'reflectedOrdinary', 'reflectedExtraordinary'})));
verifyTrue(testCase, all(isfield(interaction.Q, ...
    {'transmitted', 'reflectedOrdinary', 'reflectedExtraordinary'})));
verifyLessThan(testCase, norm(interaction.diagnostics.boundaryResidual_m), 1e-12);
verifyEqual(testCase, real(interaction.diagnostics.totalFluxRatioField), ...
    1, 'AbsTol', 1e-12);
end


function testTirDisablesTransmittedChild(testCase)
calcite = struct('nO', 1.6557, 'nE', 1.4852);
prism = createRightTriangularPrismSolid(calcite, ...
    LegY=10, LegZ=10, Width=8, OpticAxis=[1;0;0]);
options = prtDefaultOptions();
options.maxInteractions = 3;
options.minAmplitude = 1e-8;
output = polarizationRayTraceSolid( ...
    prism, [0;0;1], [0;5;-1], [1;0], options);

interactions = output.interactions( ...
    [output.interactions.caseName] == "uniaxialToIsotropic");
tirFlags = arrayfun(@(value) value.diagnostics.tir, interactions);
interaction = interactions(find(tirFlags, 1));

verifyNotEmpty(testCase, interaction);
verifyFalse(testCase, interaction.children(1).active);
verifyEqual(testCase, [interaction.children(2:3).branchType], ...
    ["reflected", "reflected"]);
verifyLessThan(testCase, norm(interaction.diagnostics.boundaryResidual_m), 1e-12);
verifyEqual(testCase, real(interaction.diagnostics.totalFluxRatioField), ...
    1, 'AbsTol', 1e-12);
end


function interaction = firstInteraction(output, caseName)
index = find([output.interactions.caseName] == caseName, 1);
interaction = output.interactions(index);
end
