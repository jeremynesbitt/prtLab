function child = makeChildTemplate(parentRay, hit, normal, mediumOut, mode)
%MAKECHILDTEMPLATE Build a child ray placeholder from a parent ray.

child = emptyRayBranch();
child.surfaceIndex = parentRay.surfaceIndex + 1;
child.mode = string(mode);
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
