% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function clearAperture = prtClearAperture(surfaceRow, defaultClearAperture)
%PRTCLEARAPERTURE Return clear aperture from SurfaceData or a default.

if nargin < 2
    defaultClearAperture = 1.0;
end

surfaceData = surfaceRow.SurfaceData{1};
if isfield(surfaceData, 'clearAperture')
    clearAperture = surfaceData.clearAperture;
elseif isfield(surfaceData, 'clear_aperture')
    clearAperture = surfaceData.clear_aperture;
else
    clearAperture = defaultClearAperture;
end

end
