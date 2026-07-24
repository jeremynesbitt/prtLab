function solid = transformPrtSolid(solid, rotation, translation)
%TRANSFORMPRTSOLID Apply a rigid transform to a prtLab solid.

arguments
    solid (1,1) struct
    rotation (3,3) double = eye(3)
    translation (3,1) double = zeros(3,1)
end

if norm(rotation.'*rotation - eye(3), 'fro') > 1e-10 || ...
        abs(det(rotation) - 1) > 1e-10
    error('transformPrtSolid:InvalidRotation', ...
        'rotation must be a proper orthonormal 3-by-3 matrix.');
end

for faceIndex = 1:numel(solid.faces)
    solid.faces(faceIndex).vertices = ...
        rotation*solid.faces(faceIndex).vertices + translation;
    solid.faces(faceIndex).point = ...
        rotation*solid.faces(faceIndex).point + translation;
    solid.faces(faceIndex).normal = rotation*solid.faces(faceIndex).normal;
end

if isfield(solid.material, 'AxisData') && ...
        isfield(solid.material.AxisData, 'opticAxis')
    solid.material.AxisData.opticAxis = ...
        rotation*solid.material.AxisData.opticAxis;
end

if isfield(solid, 'geometry') && isfield(solid.geometry, 'origin')
    solid.geometry.origin = rotation*solid.geometry.origin + translation;
end

end
