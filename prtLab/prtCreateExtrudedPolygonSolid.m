function solid = prtCreateExtrudedPolygonSolid(name, verticesYZ, faceNames, ...
        width, outside, material, lambda, lambdaUnits)
%PRTCREATEEXTRUDEDPOLYGONSOLID Build a homogeneous polygon solid extruded in x.
%
%   verticesYZ is 2-by-N with rows [y; z]. Each named face connects one
%   vertex to the next, including the final-to-first closing edge.

arguments
    name (1,1) string
    verticesYZ (2,:) double
    faceNames (:,1) string
    width (1,1) double {mustBePositive}
    outside (1,1) struct
    material (1,1) struct
    lambda (1,1) double {mustBePositive}
    lambdaUnits (1,1) string
end

numVertices = size(verticesYZ, 2);
if numVertices < 3
    error('prtCreateExtrudedPolygonSolid:TooFewVertices', ...
        'An extruded polygon requires at least three cross-section vertices.');
end
if numel(faceNames) ~= numVertices
    error('prtCreateExtrudedPolygonSolid:FaceNameCount', ...
        'Provide one face name for each polygon edge.');
end

solid = struct();
solid.type = "solid";
solid.name = name;
solid.lambda = lambda;
solid.lambdaUnits = lambdaUnits;
solid.outside = outside;
solid.material = material;
solid.faces = repmat(emptyFace(), 0, 1);

centerYZ = mean(verticesYZ, 2);
for ii = 1:numVertices
    nextIndex = mod(ii, numVertices) + 1;
    solid.faces(end+1) = makeFace( ...
        faceNames(ii), verticesYZ(:,ii), verticesYZ(:,nextIndex), ...
        width, centerYZ);
end

end

function face = emptyFace()
face = struct( ...
    'name', "", ...
    'vertices', zeros(3,0), ...
    'point', zeros(3,1), ...
    'normal', zeros(3,1), ...
    'CoatingData', struct('type', 'bare'));
end

function face = makeFace(name, p0, p1, width, centerYZ)
xMin = -width/2;
xMax = width/2;
v1 = [xMin; p0(1); p0(2)];
v2 = [xMax; p0(1); p0(2)];
v3 = [xMax; p1(1); p1(2)];
v4 = [xMin; p1(1); p1(2)];

vertices = [v1, v2, v3, v4];
normal = cross(v2 - v1, v3 - v2);
normal = normal / norm(normal);
facePoint = mean(vertices, 2);
solidCenter = [0; centerYZ(1); centerYZ(2)];
if dot(normal, facePoint - solidCenter) < 0
    vertices = fliplr(vertices);
    normal = -normal;
end

face = emptyFace();
face.name = name;
face.vertices = vertices;
face.point = mean(vertices, 2);
face.normal = normal;
end
