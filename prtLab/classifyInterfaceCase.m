% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function caseName = classifyInterfaceCase(materialIn, materialOut)
%CLASSIFYINTERFACECASE Return dispatch key for adjacent material types.

materialIn = lower(string(materialIn));
materialOut = lower(string(materialOut));

if materialIn == "isotropic" && materialOut == "isotropic"
    caseName = "isotropicToIsotropic";
elseif materialIn == "isotropic" && materialOut == "uniaxial"
    caseName = "isotropicToUniaxial";
elseif materialIn == "uniaxial" && materialOut == "isotropic"
    caseName = "uniaxialToIsotropic";
elseif materialIn == "uniaxial" && materialOut == "uniaxial"
    caseName = "uniaxialToUniaxial";
else
    caseName = "unsupported";
end

end
