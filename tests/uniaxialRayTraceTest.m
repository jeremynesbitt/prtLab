function results = uniaxialRayTraceTest(action, options)
%UNIAXIALRAYTRACETEST Regression test for isotropic-to-uniaxial tracing.
%
%   results = uniaxialRayTraceTest()
%   results = uniaxialRayTraceTest("test")
%   results = uniaxialRayTraceTest("generate")
%
%   "generate" rewrites tests/baselines/uniaxialRayTraceBaseline.csv.
%   "test" compares the current results against that baseline.

arguments
    action (1,1) string {mustBeMember(action, ["test", "generate", "table"])} = "test"
    options.Tolerance (1,1) double {mustBePositive} = 1e-10
end

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'prtLab'));

baselinePath = fullfile(repoRoot, 'tests', 'baselines', ...
    'uniaxialRayTraceBaseline.csv');

results = buildUniaxialRayTraceTable();
validatePhysicsDiagnostics(results);

switch action
    case "generate"
        baselineDir = fileparts(baselinePath);
        if ~exist(baselineDir, 'dir')
            mkdir(baselineDir);
        end
        writetable(results, baselinePath);
        fprintf('Wrote %d uniaxial ray trace baseline rows to:\n%s\n', ...
            height(results), baselinePath);

    case "test"
        if ~exist(baselinePath, 'file')
            error('uniaxialRayTraceTest:MissingBaseline', ...
                'Missing baseline file. Run uniaxialRayTraceTest("generate") first.');
        end
        baseline = readtable(baselinePath);
        compareTables(results, baseline, options.Tolerance);
        fprintf('uniaxialRayTraceTest passed: %d rows, max abs diff <= %.3g.\n', ...
            height(results), options.Tolerance);

    case "table"
        % Return only.
end

end

function results = buildUniaxialRayTraceTable()
nO = 1.666450305378000;
nE = 1.490046084791000;
lambda = 0.633;

thetaXArr = [0; 70; 35];
thetaYArr = [70; 0; 35];

axisArr = [ ...
    1, 0, 0; ...
    0, 1, 0; ...
    0, sqrt(2)/2, sqrt(2)/2; ...
    0, -sqrt(2)/2, sqrt(2)/2; ...
    0, 0, 1; ...
    1, 1, 1; ...
    1, 0, 1];

rows = [];
caseId = 1;
for axisIndex = 1:size(axisArr, 1)
    opticAxis = axisArr(axisIndex,:).';
    opticAxis = opticAxis / norm(opticAxis);

    for angleIndex = 1:numel(thetaXArr)
        thetaX = thetaXArr(angleIndex);
        thetaY = thetaYArr(angleIndex);
        kInc = calcKFromThetaXThetaY(thetaX, thetaY);

        T = makeUniaxialInterfaceSystem(lambda, nO, nE, opticAxis);
        rayOutput = polarizationRayTrace(T, kInc, [0;0;0], [1;0]);
        row = summarizeRayTrace(caseId, opticAxis, thetaX, thetaY, ...
            kInc, rayOutput);
        rows = [rows; row]; %#ok<AGROW>
        caseId = caseId + 1;
    end
end

results = struct2table(rows);
end

function T = makeUniaxialInterfaceSystem(lambda, nO, nE, opticAxis)
T = createOpticalSystem(lambda);
T.Properties.UserData.lambdaUnits = 'um';

airIndex = struct('n', 1.0);
crystalIndex = struct('nO', nO, 'nE', nE);
bareCoating = struct('type', 'bare');

T = addSurface(T, Inf, 0, "plane", struct(), ...
    "isotropic", airIndex, struct(), bareCoating);
T = addSurface(T, Inf, 0, "plane", struct(), ...
    "uniaxial", crystalIndex, struct('opticAxis', opticAxis), bareCoating);
end

function row = summarizeRayTrace(caseId, opticAxis, thetaX, thetaY, kInc, rayOutput)
interaction = rayOutput.interactions(1);
ordinary = childByMode(interaction.children, "ordinary");
extraordinary = childByMode(interaction.children, "extraordinary");
coeff = interaction.coefficients;

row = struct();
row.caseId = caseId;
row.axisX = opticAxis(1);
row.axisY = opticAxis(2);
row.axisZ = opticAxis(3);
row.thetaX_deg = thetaX;
row.thetaY_deg = thetaY;

row.kIncX = kInc(1);
row.kIncY = kInc(2);
row.kIncZ = kInc(3);

row.kOrdX = ordinary.k(1);
row.kOrdY = ordinary.k(2);
row.kOrdZ = ordinary.k(3);
row.kExtX = extraordinary.k(1);
row.kExtY = extraordinary.k(2);
row.kExtZ = extraordinary.k(3);

row.SOrdX = ordinary.S(1);
row.SOrdY = ordinary.S(2);
row.SOrdZ = ordinary.S(3);
row.SExtX = extraordinary.S(1);
row.SExtY = extraordinary.S(2);
row.SExtZ = extraordinary.S(3);

row.nOrd = ordinary.metadata.n;
row.nExt = extraordinary.metadata.n;
row.nExt_SE = extraordinary.metadata.n_SE;

row.ampOrd = ordinary.amplitude;
row.ampExt = extraordinary.amplitude;
row.fluxOrd = ordinary.flux;
row.fluxExt = extraordinary.flux;

row.a_s_to = coeff.a_s_to;
row.a_s_te = coeff.a_s_te;
row.a_s_rs = coeff.a_s_rs;
row.a_s_rp = coeff.a_s_rp;
row.a_p_to = coeff.a_p_to;
row.a_p_te = coeff.a_p_te;
row.a_p_rs = coeff.a_p_rs;
row.a_p_rp = coeff.a_p_rp;

row.boundaryResidualSNorm = norm(interaction.diagnostics.boundaryResidual_s);
row.boundaryResidualPNorm = norm(interaction.diagnostics.boundaryResidual_p);
row.kOrdSnellResidualX = ordinary.metadata.n * ordinary.k(1) - kInc(1);
row.kOrdSnellResidualY = ordinary.metadata.n * ordinary.k(2) - kInc(2);
row.kExtSnellResidualX = extraordinary.metadata.n * extraordinary.k(1) - kInc(1);
row.kExtSnellResidualY = extraordinary.metadata.n * extraordinary.k(2) - kInc(2);
row.energyDirectionSeparation = norm(ordinary.S - extraordinary.S);
row.waveNormalSeparation = norm(ordinary.k - extraordinary.k);
end

function validatePhysicsDiagnostics(results)
residualTolerance = 1e-9;

checks = { ...
    'boundaryResidualSNorm', ...
    'boundaryResidualPNorm', ...
    'kOrdSnellResidualX', ...
    'kOrdSnellResidualY', ...
    'kExtSnellResidualX', ...
    'kExtSnellResidualY'};

for ii = 1:numel(checks)
    name = checks{ii};
    maxValue = max(abs(results.(name)));
    if maxValue > residualTolerance
        error('uniaxialRayTraceTest:PhysicsDiagnosticFailed', ...
            '%s max abs value %.16g exceeds %.3g.', ...
            name, maxValue, residualTolerance);
    end
end
end

function child = childByMode(children, modeName)
modeName = string(modeName);
matches = find([children.mode] == modeName);
if numel(matches) ~= 1
    error('uniaxialRayTraceTest:UnexpectedBranchCount', ...
        'Expected exactly one "%s" branch, found %d.', modeName, numel(matches));
end
child = children(matches);
end

function compareTables(current, baseline, tolerance)
if height(current) ~= height(baseline)
    error('uniaxialRayTraceTest:RowCountMismatch', ...
        'Current row count %d differs from baseline row count %d.', ...
        height(current), height(baseline));
end

if ~isequal(current.Properties.VariableNames, baseline.Properties.VariableNames)
    error('uniaxialRayTraceTest:ColumnMismatch', ...
        'Current columns differ from baseline columns.');
end

names = current.Properties.VariableNames;
maxDiff = 0;
maxName = "";
maxRow = NaN;
for ii = 1:numel(names)
    name = names{ii};
    diffValue = abs(current.(name) - baseline.(name));
    [columnMax, rowIndex] = max(diffValue);
    if columnMax > maxDiff
        maxDiff = columnMax;
        maxName = string(name);
        maxRow = rowIndex;
    end
end

if maxDiff > tolerance
    error('uniaxialRayTraceTest:BaselineMismatch', ...
        'Baseline mismatch: max abs diff %.16g in column %s row %d exceeds tolerance %.3g.', ...
        maxDiff, maxName, maxRow, tolerance);
end
end
