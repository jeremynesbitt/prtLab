% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function ray = prtPropagateRayToHit( ...
    ray, hit, lambda, updatePosition, encodePropagationPhaseInP)
%PRTPROPAGATERAYTOHIT Accumulate segment OPL and optionally encode its phase.

arguments
    ray (1,1) struct
    hit (3,1) double
    lambda (1,1) double
    updatePosition (1,1) logical = false
    encodePropagationPhaseInP (1,1) logical = false
end

segment = hit - ray.position;
segmentOPL = rayPhaseIndex(ray) * dot(ray.k, segment);
ray.OPL = ray.OPL + segmentOPL;
ray.metadata.lastSegmentOPL = segmentOPL;

hasP = hasPropagationMatrix(ray);
if encodePropagationPhaseInP
    phase = exp(1i * 2*pi * segmentOPL / lambda);
    ray.fieldE = phase * ray.fieldE;
    ray.fieldH = phase * ray.fieldH;
    ray.E = ray.fieldE;
    ray.H = ray.fieldH;
    ray.metadata.lastSegmentPhase = phase;
    if hasP
        Popl = (ray.metadata.P_interface - ...
            ray.metadata.propagationProjector) * phase + ...
            ray.metadata.propagationProjector;
        ray.P = Popl * ray.metadata.P_beforePropagation;
        ray.metadata.P_propagated = Popl;
    end
else
    ray.metadata.lastSegmentPhase = [];
    if hasP
        ray.metadata.P_propagated = ray.metadata.P_interface;
    end
end

if updatePosition
    ray.position = hit;
end

end

function phaseIndex = rayPhaseIndex(ray)
if isfield(ray.metadata, 'phaseIndex') && ~isempty(ray.metadata.phaseIndex)
    phaseIndex = ray.metadata.phaseIndex;
elseif isfield(ray.metadata, 'n') && ~isempty(ray.metadata.n)
    phaseIndex = ray.metadata.n;
else
    phaseIndex = 1;
end
end

function tf = hasPropagationMatrix(ray)
tf = isfield(ray.metadata, 'P_interface') && ...
    isfield(ray.metadata, 'P_beforePropagation') && ...
    isfield(ray.metadata, 'propagationProjector') && ...
    ~isempty(ray.metadata.P_interface) && ...
    ~isempty(ray.metadata.P_beforePropagation) && ...
    ~isempty(ray.metadata.propagationProjector);
end
