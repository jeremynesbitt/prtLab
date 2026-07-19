function solid = createRightTriangularPrismSolid(indexData, options)
%CREATERIGHTTRIANGULARPRISMSOLID Create a homogeneous right triangular prism.
%
%   solid = createRightTriangularPrismSolid(n) creates an isotropic prism.
%
%   solid = createRightTriangularPrismSolid(struct('nO',nO,'nE',nE), ...
%       OpticAxis=[0;1;0]) creates a uniaxial crystal prism.
%
%   The right-angle vertex is options.Origin. LegY and LegZ extend from
%   that vertex in the signed global y and z directions. Width is the
%   extrusion extent in x. End caps are intentionally omitted, matching
%   createFresnelRhombSolid; Width bounds the optical faces and plotting.

arguments
    indexData
    options.LegY (1,1) double {mustBePositive} = 10
    options.LegZ (1,1) double {mustBePositive} = 10
    options.Width (1,1) double {mustBePositive} = 8
    options.Origin (3,1) double = [0;0;0]
    options.YDirection (1,1) double {mustBeMember(options.YDirection,[-1,1])} = 1
    options.ZDirection (1,1) double {mustBeMember(options.ZDirection,[-1,1])} = 1
    options.OpticAxis (3,1) double = [0;1;0]
    options.OutsideIndex (1,1) double {mustBePositive} = 1.0
    options.lambda (1,1) double {mustBePositive} = 0.633
    options.lambdaUnits (1,1) string = "um"
    options.Name (1,1) string = "Right triangular prism"
end

outside = isotropicMedium(options.OutsideIndex);
material = materialMedium(indexData, options.OpticAxis);

y0 = options.Origin(2);
z0 = options.Origin(3);
yLeg = options.YDirection * options.LegY;
zLeg = options.ZDirection * options.LegZ;

% A is the right-angle vertex. Edge A-B is normal to z, B-C is the
% hypotenuse, and C-A is normal to y.
A = [y0; z0];
B = A + [yLeg; 0];
C = A + [0; zLeg];
verticesYZ = [A, B, C];
faceNames = ["z leg"; "hypotenuse"; "y leg"];

solid = prtCreateExtrudedPolygonSolid(options.Name, verticesYZ, faceNames, ...
    options.Width, outside, material, options.lambda, options.lambdaUnits);
solid.geometry = struct( ...
    'shape', "rightTriangle", ...
    'origin', options.Origin, ...
    'legY', options.LegY, ...
    'legZ', options.LegZ, ...
    'yDirection', options.YDirection, ...
    'zDirection', options.ZDirection, ...
    'width', options.Width);

end

function medium = materialMedium(indexData, opticAxis)
if isnumeric(indexData) && isscalar(indexData)
    medium = isotropicMedium(indexData);
    return;
end
if ~isstruct(indexData) || ~isscalar(indexData)
    error('createRightTriangularPrismSolid:InvalidIndexData', ...
        'indexData must be a scalar refractive index or a scalar struct.');
end

if isfield(indexData, 'n')
    medium = struct( ...
        'MaterialType', "isotropic", ...
        'IndexData', indexData, ...
        'AxisData', struct());
elseif isfield(indexData, 'nO') && isfield(indexData, 'nE')
    if norm(opticAxis) == 0
        error('createRightTriangularPrismSolid:ZeroOpticAxis', ...
            'OpticAxis must be nonzero for a uniaxial material.');
    end
    medium = struct( ...
        'MaterialType', "uniaxial", ...
        'IndexData', indexData, ...
        'AxisData', struct('opticAxis', opticAxis/norm(opticAxis)));
else
    error('createRightTriangularPrismSolid:MissingIndices', ...
        'indexData must contain n, or both nO and nE.');
end
end

function medium = isotropicMedium(n)
medium = struct( ...
    'MaterialType', "isotropic", ...
    'IndexData', struct('n', n), ...
    'AxisData', struct());
end
