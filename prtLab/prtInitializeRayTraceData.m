% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function rayTraceData = prtInitializeRayTraceData(system, options)
%PRTINITIALIZERAYTRACEDATA Create the common ray-trace output container.

rayTraceData = struct();
rayTraceData.system = system;
rayTraceData.options = options;
rayTraceData.rays = repmat(emptyRayBranch(), 0, 1);
rayTraceData.interactions = repmat(emptyInteractionRecord(), 0, 1);
rayTraceData.finalRayIds = [];

end
