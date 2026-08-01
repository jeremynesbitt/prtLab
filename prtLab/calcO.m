% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function Oout = calcO(s,p,k)

assert((size(s,1) == 3) && ...
       (size(s,2) == 1) && ...
       (all(size(s) == size(p))) && ...
       (all(size(s) == size(k))),'Error!  vectors must all be 3x1');

% All column vectors
Oout = [s, p, k];


end