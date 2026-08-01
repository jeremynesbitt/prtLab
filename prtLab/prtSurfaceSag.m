% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function zSag = prtSurfaceSag(x, y, surfaceRow)
%PRTSURFACESAG Evaluate sag for supported prtLab surface rows.

surfaceType = string(surfaceRow.SurfaceType);
switch surfaceType
    case "plane"
        zSag = 0;

    case "odd_asphere"
        surfaceData = surfaceRow.SurfaceData{1};
        [kappa, A] = oddAsphereParameters(surfaceData);
        zSag = oddAsphereSag(x, y, surfaceRow.Radius, kappa, A);

    otherwise
        error('prtSurfaceSag:UnsupportedSurface', ...
            'Surface type "%s" is not implemented.', surfaceType);
end

end

function [kappa, A] = oddAsphereParameters(surfaceData)
if isfield(surfaceData, 'kappa')
    kappa = surfaceData.kappa;
elseif isfield(surfaceData, 'data')
    kappa = surfaceData.data(1);
else
    kappa = 0;
end

if isfield(surfaceData, 'A')
    A = surfaceData.A(:).';
elseif isfield(surfaceData, 'coefficients')
    A = surfaceData.coefficients(:).';
elseif isfield(surfaceData, 'data')
    A = surfaceData.data(2:end);
else
    A = zeros(1,8);
end

if numel(A) < 8
    A = [A, zeros(1, 8-numel(A))];
elseif numel(A) > 8
    A = A(1:8);
end
end

function zSag = oddAsphereSag(x, y, radius, kappa, A)
r2 = x^2 + y^2;
r = sqrt(r2);
if isinf(radius)
    zBase = 0;
else
    c = 1 / radius;
    arg = max(1 - kappa*c^2*r2, 0);
    zBase = c*r2 / (1 + sqrt(arg));
end

zPoly = 0;
for power = 3:10
    zPoly = zPoly + A(power-2) * r^power;
end
zSag = zBase + zPoly;
end
