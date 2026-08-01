% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [hit, normal, interceptData] = intersectSurface(surfaceRow, ray)
%INTERSECTSURFACE Compute ray intercept and normal for one surface row.
%
%   Surface locations are derived sequentially from the row Thickness:
%   target vertex z = ray.position(3) + Thickness. Odd-asphere surfaces use
%   the Chipman Eq. 10.65 convention, with odd and even powers A3..A10.

surfaceType = string(surfaceRow.SurfaceType);
thickness = surfaceRow.Thickness;
surfaceData = surfaceRow.SurfaceData{1};
previousVertexZ = currentVertexZ(ray);
targetVertexZ = previousVertexZ + thickness;

switch surfaceType
    case "plane"
        targetZ = targetVertexZ;
        step = (targetZ - ray.position(3)) / ray.S(3);
        hit = ray.position + step * ray.S;
        normal = [0; 0; 1];

        interceptData = struct( ...
            'step', step, ...
            'targetZ', targetZ, ...
            'surfaceType', surfaceType);

    case "odd_asphere"
        radius = surfaceRow.Radius;
        vertexZ = targetVertexZ;
        [kappa, A] = oddAsphereParameters(surfaceData);
        [hit, step, iterations, residual] = oddAsphereIntercept( ...
            ray.position, ray.S, vertexZ, radius, kappa, A);
        normal = oddAsphereNormal(hit(1), hit(2), radius, kappa, A);

        interceptData = struct( ...
            'step', step, ...
            'targetZ', vertexZ, ...
            'surfaceType', surfaceType, ...
            'kappa', kappa, ...
            'A', A, ...
            'iterations', iterations, ...
            'residual', residual);

    otherwise
        error('intersectSurface:UnsupportedSurface', ...
            'Surface type "%s" is not implemented yet.', surfaceType);
end

normal = normal / norm(normal);
if dot(normal, ray.S) < 0
    normal = -normal;
end

end

function vertexZ = currentVertexZ(ray)
if isfield(ray.metadata, 'currentVertexZ')
    vertexZ = ray.metadata.currentVertexZ;
else
    vertexZ = 0;
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

function [hit, step, iterations, residual] = oddAsphereIntercept(r0, k, vertexZ, radius, kappa, A)
if abs(k(3)) > 1e-12
    step = (vertexZ - r0(3)) / k(3);
else
    step = 1e-3;
end

ep = 1e-7;
residual = NaN;
for iterations = 1:300
    x = r0(1) + step*k(1);
    y = r0(2) + step*k(2);
    z = r0(3) + step*k(3);
    residual = z - vertexZ - oddAsphereSag(x, y, radius, kappa, A);
    dsdx = (oddAsphereSag(x+ep, y, radius, kappa, A) - ...
        oddAsphereSag(x-ep, y, radius, kappa, A)) / (2*ep);
    dsdy = (oddAsphereSag(x, y+ep, radius, kappa, A) - ...
        oddAsphereSag(x, y-ep, radius, kappa, A)) / (2*ep);
    dfdt = k(3) - dsdx*k(1) - dsdy*k(2);
    if abs(dfdt) < 1e-15
        break;
    end
    deltaStep = -residual / dfdt;
    step = step + deltaStep;
    if abs(deltaStep) < 1e-12
        break;
    end
end

hit = r0 + step*k;
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

function normal = oddAsphereNormal(x, y, radius, kappa, A)
ep = 1e-7;
dsdx = (oddAsphereSag(x+ep, y, radius, kappa, A) - ...
    oddAsphereSag(x-ep, y, radius, kappa, A)) / (2*ep);
dsdy = (oddAsphereSag(x, y+ep, radius, kappa, A) - ...
    oddAsphereSag(x, y-ep, radius, kappa, A)) / (2*ep);
normal = [-dsdx; -dsdy; 1];
normal = normal / norm(normal);
end
