% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function tf = isTotalInternalReflection(k, normal, n1, n2)
cosThetaI = dot(k, normal);
sin2ThetaI = max(0, 1 - cosThetaI^2);
tf = isreal(n1) && isreal(n2) && n1 > n2 && (n1/n2)^2 * sin2ThetaI > 1;
end
