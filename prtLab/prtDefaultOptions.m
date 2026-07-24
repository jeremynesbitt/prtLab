function options = prtDefaultOptions(userOptions)
%PRTDEFAULTOPTIONS Default options for prtLab ray tracing.

if nargin < 1
    userOptions = struct();
end

options = struct();
options.maxBranches = 256;
options.minFlux = 0;
options.minRelativeFlux = 0;
options.minAmplitude = 0;
options.keepDiagnostics = true;
options.surfaceTolerance = 1e-10;
options.maxInterceptIterations = 25;
options.dispatchUnimplemented = "error";

userFields = fieldnames(userOptions);
for ii = 1:numel(userFields)
    options.(userFields{ii}) = userOptions.(userFields{ii});
end

end
