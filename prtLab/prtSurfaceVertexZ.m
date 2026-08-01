% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function vertexZ = prtSurfaceVertexZ(T)
%PRTSURFACEVERTEXZ Cumulative vertex z positions from a prtLab table.

vertexZ = cumsum(T.Thickness(:)).';

end
