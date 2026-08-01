% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [thetaX, thetaY] = getThetaXThetaYFromK(kinc)
%GETTHETAXTHETAYFROMK Invert the PLAOS thetaX/thetaY ray convention.
%
%   [thetaX, thetaY] = getThetaXThetaYFromK(kinc) returns angles in radians
%   for a direction vector parameterized as
%
%       k = [tan(thetaX); tan(thetaY); 1] / norm(...)
%
%   kinc may be any nonzero scalar multiple of this direction, but kinc(3)
%   must be nonzero.

arguments
    kinc (3,1) double
end

assert(norm(kinc) > 0, ...
    'getThetaXThetaYFromK:ZeroVector', ...
    'kinc must be nonzero.');

assert(abs(kinc(3)) > eps(norm(kinc)), ...
    'getThetaXThetaYFromK:GrazingRay', ...
    'Cannot recover finite thetaX/thetaY when kinc(3) is zero.');

thetaX = atan(kinc(1) / kinc(3));
thetaY = atan(kinc(2) / kinc(3));

end
