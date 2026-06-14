function vertexZ = prtSurfaceVertexZ(T)
%PRTSURFACEVERTEXZ Cumulative vertex z positions from a prtLab table.

vertexZ = cumsum(T.Thickness(:)).';

end
