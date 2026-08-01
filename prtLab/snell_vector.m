% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function k_t = snell_vector(k_i, ada, n1, n2)
% Vector form of Snell's law (Eq. 10.20 from Chipman)
% k_i: incident unit vector, 
% ada: surface normal into transmitted medium (unit)
% n1: incident index 
% n2: transmitted index
% Rev 2 - support complex k (tir use case)

cos_i = dot(k_i, ada);
nr = n1/n2;
disc = 1 - nr^2*(1 - cos_i^2);

if isreal(disc) && disc < 0
    root = 1i*sqrt(-disc);
else
    root = sqrt(disc);
end

k_t = nr*k_i - (nr*cos_i - root)*ada; % same as Chipman 10.20 

k_t = prtNorm(k_t);

end
