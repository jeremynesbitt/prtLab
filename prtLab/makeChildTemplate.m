% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function child = makeChildTemplate(parentRay, hit, normal, mediumOut, mode, branchType)
%MAKECHILDTEMPLATE Build a child ray placeholder from a parent ray.

arguments
    parentRay (1,1) struct
    hit (3,1) double
    normal (3,1) double
    mediumOut
    mode (1,1) string
    branchType (1,1) string {mustBeMember(branchType, ...
        ["transmitted", "reflected"])}
end

child = emptyRayBranch();
child.surfaceIndex = parentRay.surfaceIndex + 1;
child.mode = mode;
child.branchType = branchType;
child.mediumType = string(mediumOut.MaterialType);
child.position = hit;
child.k = parentRay.k;
child.S = parentRay.S;
child.modeE = parentRay.modeE;
child.modeH = parentRay.modeH;
child.fieldE = parentRay.fieldE;
child.fieldH = parentRay.fieldH;
child.E = parentRay.E;
child.H = parentRay.H;
child.P = parentRay.P;
child.Q = parentRay.Q;
child.O = parentRay.O;
child.localBasis = parentRay.localBasis;
child.amplitude = parentRay.amplitude;
child.flux = parentRay.flux;
child.OPL = parentRay.OPL;
child.active = true;
child.metadata = struct('normal', normal);

end
