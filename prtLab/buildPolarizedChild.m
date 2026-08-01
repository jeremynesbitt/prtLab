% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function child = buildPolarizedChild(parentRay, hit, normal, medium, options)
%BUILDPOLARIZEDCHILD Build a complete ray child from modal branch data.

arguments
    parentRay (1,1) struct
    hit (3,1) double
    normal (3,1) double
    medium
    options.Mode (1,1) string
    options.BranchType (1,1) string {mustBeMember(options.BranchType, ...
        ["transmitted", "reflected"])}
    options.k (3,1) double
    options.S (3,1) double
    options.P (3,3) double
    options.Q (3,3) double
    options.O (3,3) double
    options.LocalBasis (1,1) struct
    options.Index (1,1) double
    options.ModeE double = []
    options.ModeH double = []
    options.FluxNormal (3,1) double = normal
    options.Active (1,1) logical = true
    options.PhaseIndex double = []
    options.PropagationProjector double = []
    options.Metadata (1,1) struct = struct()
end

child = makeChildTemplate(parentRay, hit, normal, medium, ...
    options.Mode, options.BranchType);
child.k = options.k;
child.S = options.S;
child.P = options.P * parentRay.P;
child.Q = options.Q * parentRay.Q;
child.O = options.O;
child.localBasis = options.LocalBasis;
child.active = options.Active;

child.fieldE = options.P * parentRay.fieldE;
child.modeE = resolveModeE(options.ModeE, child.fieldE);
child.modeH = resolveModeH(options.ModeH, child.modeE, options.Index, options.k);
modalScale = prtModalFieldScale(child.modeE, child.fieldE);
child.fieldH = modalScale * child.modeH;
child.E = child.fieldE;
child.H = child.fieldH;
child.amplitude = norm(child.fieldE);
child.flux = real(dot(cross(child.fieldE, conj(child.fieldH)), ...
    options.FluxNormal));

metadata = options.Metadata;
phaseIndex = options.PhaseIndex;
if isempty(phaseIndex)
    phaseIndex = options.Index;
end
metadata.normal = normal;
metadata.n = options.Index;
metadata.phaseIndex = phaseIndex;
metadata.modalScale = modalScale;
metadata.P_interface = options.P;
metadata.P_beforePropagation = parentRay.P;
metadata.IoutField = child.flux;
if ~isempty(options.PropagationProjector)
    metadata.propagationProjector = options.PropagationProjector;
end
child.metadata = metadata;

end


function modeE = resolveModeE(requestedModeE, fieldE)
if ~isempty(requestedModeE)
    modeE = requestedModeE(:);
elseif norm(fieldE) > 0
    modeE = prtNorm(fieldE);
else
    modeE = zeros(3,1);
end
end


function modeH = resolveModeH(requestedModeH, modeE, index, k)
if ~isempty(requestedModeH)
    modeH = requestedModeH(:);
elseif norm(modeE) > 0
    modeH = index * makeK(k) * modeE;
else
    modeH = zeros(3,1);
end
end
