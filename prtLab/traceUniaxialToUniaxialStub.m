% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function interaction = traceUniaxialToUniaxialStub(interaction, mediumIn, mediumOut, ray, hit, normal, options)
%TRACEUNIAXIALTOUNIAXIALSTUB Placeholder for uniaxial-to-uniaxial tracing.

if options.dispatchUnimplemented == "error"
    error('traceUniaxialToUniaxialStub:NotImplemented', ...
        'uniaxial-to-uniaxial tracing is not implemented yet.');
end

childO = makeChildTemplate(ray, hit, normal, mediumOut, ...
    "ordinary", "transmitted");
childE = makeChildTemplate(ray, hit, normal, mediumOut, ...
    "extraordinary", "transmitted");
interaction.children = [childO; childE];

end
