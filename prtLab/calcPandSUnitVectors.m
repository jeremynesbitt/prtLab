% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [p, s] = calcPandSUnitVectors(k, ada, fallbackS)
%CALCPANDSUNITVECTORS Build p/s unit vectors for a ray and interface normal.
%
%   [p, s] = calcPandSUnitVectors(k, ada)
%   [p, s] = calcPandSUnitVectors(k, ada, fallbackS)
%
%   s is parallel to cross(k, ada). If k and ada are parallel, fallbackS is
%   projected into the transverse plane and used instead. Without fallbackS,
%   a stable global-axis fallback is chosen.

arguments
    k (3,1) double
    ada (3,1) double
    fallbackS (3,1) double = defaultFallbackS(k)
end

k = prtNorm(k);
ada = prtNorm(ada);

s = cross(k, ada);
if norm(s) < 100*eps(max(1, norm(k)*norm(ada)))
    s = fallbackS - dot(fallbackS, k)*k;
    if norm(s) < 100*eps(max(1, norm(fallbackS)))
        s = defaultFallbackS(k);
        s = s - dot(s, k)*k;
    end
end
s = prtNorm(s);

p = cross(k, s);
p = prtNorm(p);
end

function s = defaultFallbackS(k)
if abs(dot(k(:)/norm(k), [1;0;0])) < 0.9
    s = [1;0;0];
else
    s = [0;1;0];
end
end
