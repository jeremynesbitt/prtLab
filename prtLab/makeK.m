% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function K = makeK(k_inc)

K = [0, -k_inc(3), k_inc(2); ...
        k_inc(3), 0, -k_inc(1) ; ...
        -k_inc(2), k_inc(1), 0];

end