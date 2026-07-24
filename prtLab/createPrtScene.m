function scene = createPrtScene(lambda, outsideMedium, options)
%CREATEPRTSCENE Create a scene containing multiple homogeneous solids.

arguments
    lambda (1,1) double {mustBePositive}
    outsideMedium = 1.0
    options.lambdaUnits (1,1) string = "um"
    options.Name (1,1) string = "prtLab scene"
end

scene = struct();
scene.type = "scene";
scene.name = options.Name;
scene.lambda = lambda;
scene.lambdaUnits = options.lambdaUnits;
scene.outside = normalizeMedium(outsideMedium);
scene.solids = cell(0,1);

end

function medium = normalizeMedium(value)
if isnumeric(value) && isscalar(value) && isreal(value) && value > 0
    medium = struct( ...
        'MaterialType', "isotropic", ...
        'IndexData', struct('n', value), ...
        'AxisData', struct());
    return;
end

if ~isstruct(value) || ~isscalar(value) || ...
        ~all(isfield(value, {'MaterialType', 'IndexData', 'AxisData'}))
    error('createPrtScene:InvalidOutsideMedium', ...
        ['outsideMedium must be a positive scalar index or a medium struct ', ...
         'with MaterialType, IndexData, and AxisData.']);
end
medium = value;
end
