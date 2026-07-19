%OPTIMIZEACHROMATICWAVEPLATETHICKNESS
% Find quartz/MgF2 thicknesses for a broadband quarter-wave retarder.
%
% Length units are um. Wavelengths are specified in nm for readability and
% converted to um for prtLab.

addpath(genpath('..'));

wlNm = linspace(400, 700, 17);
targetRetardanceWaves = 0.25;

T0 = exampleAchromaticWavePlateSystem();
tQuartz0 = T0.Thickness(2);
tGap = T0.Thickness(3);
tMgF20 = T0.Thickness(4);
opticAxisQuartz = T0.AxisData{2}.opticAxis;
opticAxisMgF2 = T0.AxisData{4}.opticAxis;

kIn = [0; 0; 1];
xIn = [0; 0; 0];
Ein = [1; 0];
coord = struct( ...
    'type', "doublePole", ...
    'a_loc', kIn, ...
    'x_o', [1; 0; 0]);

if abs(dot(opticAxisQuartz, opticAxisMgF2)) > 1e-10
    warning('Waveplate optic axes are not crossed. abs(dot(axis1,axis2)) = %.3g', ...
        abs(dot(opticAxisQuartz, opticAxisMgF2)));
end

x0 = [tQuartz0; tMgF20];
lb = [1; 1];
ub = [2000; 2000];

objective = @(x) retardanceResiduals( ...
    x, wlNm, targetRetardanceWaves, tGap, ...
    opticAxisQuartz, opticAxisMgF2, kIn, xIn, Ein, coord);

opts = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'FunctionTolerance', 1e-12, ...
    'StepTolerance', 1e-12);

[xOpt, resnorm, residuals, exitflag, output] = lsqnonlin(objective, x0, lb, ub, opts);

[retardanceWaves, scalarEstimateWaves] = evaluateRetardanceSpectrum( ...
    xOpt, wlNm, tGap, opticAxisQuartz, opticAxisMgF2, kIn, xIn, Ein, coord);

resultTable = table( ...
    wlNm(:), ...
    retardanceWaves(:), ...
    retardanceWaves(:) - targetRetardanceWaves, ...
    scalarEstimateWaves(:), ...
    'VariableNames', {'Wavelength_nm', 'PRT_Retardance_waves', 'Error_waves', 'ScalarEstimate_waves'});

fprintf('\nOptimized thicknesses:\n');
fprintf('  quartz: %.9g um\n', xOpt(1));
fprintf('  MgF2:   %.9g um\n', xOpt(2));
fprintf('  gap:    %.9g um (held fixed)\n', tGap);
fprintf('resnorm: %.6g, exitflag: %d\n\n', resnorm, exitflag);
disp(resultTable);
disp(output);

figure;
plot(wlNm, retardanceWaves, 'o-', 'LineWidth', 1.5);
hold on;
yline(targetRetardanceWaves, 'k--', 'LineWidth', 1.0);
plot(wlNm, scalarEstimateWaves, 's:', 'LineWidth', 1.2);
grid on;
xlabel('Wavelength [nm]');
ylabel('Retardance [waves]');
legend('prtLab Jones eigen-retardance', 'target', 'scalar crossed-axis estimate', ...
    'Location', 'best');
title('Optimized quartz/MgF2 achromatic waveplate');

function residuals = retardanceResiduals(x, wlNm, targetRetardanceWaves, tGap, ...
    opticAxisQuartz, opticAxisMgF2, kIn, xIn, Ein, coord)

retardanceWaves = evaluateRetardanceSpectrum( ...
    x, wlNm, tGap, opticAxisQuartz, opticAxisMgF2, kIn, xIn, Ein, coord);

residuals = retardanceWaves(:) - targetRetardanceWaves;
end

function [retardanceWaves, scalarEstimateWaves] = evaluateRetardanceSpectrum( ...
    x, wlNm, tGap, opticAxisQuartz, opticAxisMgF2, kIn, xIn, Ein, coord)

tQuartz = x(1);
tMgF2 = x(2);
retardanceWaves = nan(size(wlNm));
scalarEstimateWaves = nan(size(wlNm));

for ii = 1:numel(wlNm)
    lambdaUm = wlNm(ii) / 1000;
    [nOQuartz, nEQuartz] = quartzOpticalConstants(lambdaUm);
    [nOMgF2, nEMgF2] = mgF2OpticalConstants(lambdaUm);

    T = modAchromaticWavePlateSystem( ...
        lambdaUm, ...
        nOQuartz, nEQuartz, ...
        nOMgF2, nEMgF2, ...
        tQuartz, tGap, tMgF2, ...
        opticAxisQuartz, opticAxisMgF2);

    rayOutput = polarizationRayTrace(T, kIn, xIn, Ein);
    Ptot = coherentFinalP(rayOutput);
    J = transformPtoJones(Ptot, kIn, kIn, coord);

    retardanceWaves(ii) = jonesRetardanceWaves(J);
    scalarEstimateWaves(ii) = scalarCrossedAxisRetardance( ...
        lambdaUm, nOQuartz, nEQuartz, nOMgF2, nEMgF2, tQuartz, tMgF2);
end
end

function Ptot = coherentFinalP(rayOutput)
Ptot = zeros(3,3);
for kk = 1:numel(rayOutput.finalRayIds)
    Ptot = Ptot + rayOutput.rays(rayOutput.finalRayIds(kk)).P;
end
end

function retardanceWaves = jonesRetardanceWaves(J)
ev = eig(J);
if any(abs(ev) < 1e-14)
    retardanceWaves = NaN;
    return;
end

retardanceWaves = angle(ev(2) / ev(1)) / (2*pi);
retardanceWaves = mod(retardanceWaves, 1);
retardanceWaves = min(retardanceWaves, 1 - retardanceWaves);
end

function retardanceWaves = scalarCrossedAxisRetardance( ...
    lambdaUm, nOQuartz, nEQuartz, nOMgF2, nEMgF2, tQuartz, tMgF2)

deltaQuartz = (nEQuartz - nOQuartz) * tQuartz;
deltaMgF2 = (nEMgF2 - nOMgF2) * tMgF2;
retardanceWaves = mod((deltaMgF2 - deltaQuartz) / lambdaUm, 1);
retardanceWaves = min(retardanceWaves, 1 - retardanceWaves);
end
