% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function ray = prtMakeInitialRay(kInc, xInc, Ein, mediumType, indexData, metadata)
%PRTMAKEINITIALRAY Create a common input ray for a prtLab trace.

arguments
    kInc (3,1) double
    xInc (3,1) double
    Ein double
    mediumType (1,1) string
    indexData (1,1) struct
    metadata (1,1) struct = struct()
end

kInc = kInc / norm(kInc);
[fieldE, localBasis] = initializeInputElectricField(Ein, kInc);
nInc = incidentIndex(indexData);

ray = emptyRayBranch();
ray.id = 1;
ray.parentId = 0;
ray.surfaceIndex = 0;
if mediumType == "isotropic"
    ray.mode = "isotropic";
else
    ray.mode = "input";
end
ray.branchType = "input";
ray.mediumType = mediumType;
ray.position = xInc;
ray.k = kInc;
ray.S = kInc;
ray.modeE = fieldE / norm(fieldE);
ray.modeH = nInc * makeK(kInc) * ray.modeE;
ray.fieldE = fieldE;
ray.fieldH = nInc * makeK(kInc) * fieldE;
ray.E = ray.fieldE;
ray.H = ray.fieldH;
ray.P = eye(3);
ray.Q = eye(3);
ray.O = eye(3);
ray.localBasis = localBasis;
ray.amplitude = norm(fieldE);
ray.flux = real(dot(cross(fieldE, conj(ray.fieldH)), kInc));
ray.OPL = 0;
ray.active = true;
ray.history = 1;
metadata.n = nInc;
ray.metadata = metadata;

end

function [fieldE, localBasis] = initializeInputElectricField(Ein, kInc)
if numel(Ein) == 2
    Ein = Ein(:);
    [xDp, yDp] = doublePoleBasisVectors(kInc, [0;0;1], [1;0;0]);
    xDp = xDp(:);
    yDp = yDp(:);
    fieldE = Ein(1)*xDp + Ein(2)*yDp;
    localBasis = struct( ...
        'x', xDp, ...
        'y', yDp, ...
        'basisDirection', kInc, ...
        'inputConvention', "doublePoleJones");
elseif numel(Ein) == 3
    fieldE = Ein(:);
    longitudinal = dot(fieldE, kInc);
    if abs(longitudinal) > 100*eps(max(1, norm(fieldE)))
        fieldE = fieldE - longitudinal*kInc;
    end
    fieldE = fieldE(:);
    localBasis = struct( ...
        'x', [], ...
        'y', [], ...
        'basisDirection', kInc, ...
        'inputConvention', "global3D");
else
    error('prtMakeInitialRay:InvalidInputField', ...
        'Ein must be either a 2x1 Jones vector or a 3x1 electric field.');
end

if norm(fieldE) == 0
    error('prtMakeInitialRay:ZeroInputField', ...
        'Ein must have nonzero electric-field amplitude.');
end
end

function nInc = incidentIndex(indexData)
if isfield(indexData, 'n')
    nInc = indexData.n;
elseif isfield(indexData, 'nO')
    nInc = indexData.nO;
else
    error('prtMakeInitialRay:MissingIndex', ...
        'The input medium IndexData must contain n or nO.');
end
end
