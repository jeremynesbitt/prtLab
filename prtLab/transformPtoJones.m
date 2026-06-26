function Jout = transformPtoJones(P, ki, kj, coordIn, coordOut)
%TRANSFORMPTOJONES Convert a 3D P matrix to a 2D Jones matrix.
%
%   Jout = transformPtoJones(P, ki, kj, coordIn)
%   Jout = transformPtoJones(P, ki, kj, coordIn, coordOut)
%
%   P maps the 3D incident field to the 3D exiting field. ki and kj are the
%   incident and exiting propagation directions. coordIn and coordOut define
%   the entrance and exit local coordinate systems. If coordOut is omitted,
%   coordIn is used for both pupils.
%
%   Supported coordinate structs:
%
%       coord.type = "doublePole";
%       coord.a_loc = [0;0;1];
%       coord.x_o = [1;0;0];
%
%       coord.type = "dipole";
%       coord.a_loc = [0;0;1];
%
%       coord.type = "sp";
%       coord.normal = surfaceNormal;

arguments
    P (3,3) double
    ki (3,1) double
    kj (3,1) double
    coordIn (1,1) struct
    coordOut (1,1) struct = coordIn
end

ki = ki / norm(ki);
kj = kj / norm(kj);

[xi, yi] = localBasis(ki, coordIn);
[xj, yj] = localBasis(kj, coordOut);

Uin = [xi, yi, ki];
Uout = [xj, yj, kj];

J3 = Uout \ (P * Uin);
Jout = J3(1:2, 1:2);

end

function [xHat, yHat] = localBasis(k, coord)
if ~isfield(coord, 'type')
    error('transformPtoJones:MissingCoordinateType', ...
        'Coordinate struct must contain a type field.');
end

coordType = string(coord.type);
switch lower(coordType)
    case "doublepole"
        requireVectorField(coord, 'a_loc', coordType);
        requireVectorField(coord, 'x_o', coordType);
        [xRow, yRow] = doublePoleBasisVectors(k, coord.a_loc(:), coord.x_o(:));
        xHat = xRow(:);
        yHat = yRow(:);

    case "dipole"
        requireVectorField(coord, 'a_loc', coordType);
        [xRow, yRow] = dipoleBasisVectors(k(:).', coord.a_loc(:).');
        xHat = xRow(:);
        yHat = yRow(:);

    case "sp"
        requireVectorField(coord, 'normal', coordType);
        [xHat, yHat] = spBasis(k, coord.normal(:));

    otherwise
        error('transformPtoJones:UnknownCoordinateType', ...
            'Unknown coord.type "%s". Use "doublePole", "dipole", or "sp".', coordType);
end

end

function requireVectorField(coord, fieldName, coordType)
if ~isfield(coord, fieldName)
    error('transformPtoJones:MissingCoordinateField', ...
        'Coordinate type "%s" requires field "%s".', coordType, fieldName);
end

v = coord.(fieldName);
if isempty(v) || numel(v) ~= 3 || norm(v(:)) == 0
    error('transformPtoJones:InvalidCoordinateField', ...
        'Field "%s" must be a nonzero 3-vector for coordinate type "%s".', ...
        fieldName, coordType);
end
end

function [sHat, pHat] = spBasis(k, normal)
normal = normal / norm(normal);
sHat = cross(k, normal);
if norm(sHat) < 1e-12
    error('transformPtoJones:SingularSPBasis', ...
        's/p coordinates are singular because k is parallel to the surface normal.');
end
sHat = sHat / norm(sHat);
pHat = cross(k, sHat);
pHat = pHat / norm(pHat);
end
