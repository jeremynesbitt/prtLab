function rayTraceData = prtInitializeRayTraceData(system, options)
%PRTINITIALIZERAYTRACEDATA Create the common ray-trace output container.

rayTraceData = struct();
rayTraceData.system = system;
rayTraceData.options = options;
rayTraceData.rays = repmat(emptyRayBranch(), 0, 1);
rayTraceData.interactions = repmat(emptyInteractionRecord(), 0, 1);
rayTraceData.finalRayIds = [];

end
