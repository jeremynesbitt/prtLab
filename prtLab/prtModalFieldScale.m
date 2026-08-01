% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function scale = prtModalFieldScale(modeE, fieldE)
%PRTMODALFIELDSCALE Least-squares scalar mapping a modal E field to fieldE.

denom = modeE' * modeE;
if abs(denom) < eps
    scale = 0;
else
    scale = (modeE' * fieldE) / denom;
end

end
