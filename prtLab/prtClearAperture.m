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
