% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [ordinary, extraordinary] = uniaxialModesFromTangentialQ( ...
    qTangential, normal, nO, nE, opticAxis, directionSign)
%UNIAXIALMODESFROMTANGENTIALQ Solve both modes with conserved tangential q.

if nargin < 6
    directionSign = 1;
end

ada = normal(:) / norm(normal);
axis = opticAxis(:) / norm(opticAxis);
qTangential = qTangential(:);
epsilon = nO^2*eye(3) + (nE^2 - nO^2)*(axis*axis.');

ordinaryRoot = sqrt(complex(nO^2 - dot(qTangential, qTangential)));
ordinary = selectMode([ ...
    modeFromQ(qTangential + ordinaryRoot*ada, epsilon); ...
    modeFromQ(qTangential - ordinaryRoot*ada, epsilon)], ...
    ada, directionSign);

axisNormal = dot(axis, ada);
axisTangential = dot(axis, qTangential);
A = axisNormal^2/nO^2 + (1 - axisNormal^2)/nE^2;
B = 2*axisTangential*axisNormal*(1/nO^2 - 1/nE^2);
C = axisTangential^2/nO^2 + ...
    (dot(qTangential, qTangential) - axisTangential^2)/nE^2 - 1;
normalRoots = roots([A, B, C]);
extraordinary = selectMode([ ...
    modeFromQ(qTangential + normalRoots(1)*ada, epsilon); ...
    modeFromQ(qTangential + normalRoots(2)*ada, epsilon)], ...
    ada, directionSign);
end


function mode = modeFromQ(q, epsilon)
M = epsilon + makeK(q)*makeK(q);
[~,~,V] = svd(M);
E = canonicalizeMode(V(:,end));
H = cross(q, E);
poynting = real(cross(E, conj(H)));
S = poynting / norm(poynting);
phaseIndex = sqrt(sum(q.*q));
mode = struct( ...
    'q', q, ...
    'k', q/phaseIndex, ...
    'E', E, ...
    'H', H, ...
    'S', S, ...
    'phaseIndex', phaseIndex);
end


function mode = selectMode(candidates, normal, directionSign)
scores = arrayfun(@(value) ...
    directionSign*real(dot(value.S, normal)), candidates);
valid = find(scores > 0);
if isempty(valid)
    error('uniaxialModesFromTangentialQ:NoPropagatingMode', ...
        'No uniaxial mode propagates in the requested direction.');
end
[~, localIndex] = max(scores(valid));
mode = candidates(valid(localIndex));
end


function value = canonicalizeMode(value)
value = value / norm(value);
[~, pivot] = max(abs(value));
if abs(value(pivot)) > 0
    value = value * exp(-1i*angle(value(pivot)));
end
if real(value(pivot)) < 0
    value = -value;
end
if max(abs(imag(value))) < 100*eps(max(1, norm(value)))
    value = real(value);
end
end
