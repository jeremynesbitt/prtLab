% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [nO, nE] = mgF2OpticalConstants(lambda)
%lambda is in um
no2 = 1 ...
    + (0.48755108*lambda.^2)./(lambda.^2 - 0.04338408^2) ...
    + (0.39875031*lambda.^2)./(lambda.^2 - 0.09461442^2) ...
    + (2.3120353*lambda.^2)./(lambda.^2 - 23.793604^2);

nO = sqrt(no2);

ne2 = 1 ...
    + (0.41344023*lambda.^2)./(lambda.^2 - 0.03684262^2) ...
    + (0.50497499*lambda.^2)./(lambda.^2 - 0.09076162^2) ...
    + (2.4904862*lambda.^2)./(lambda.^2 - 23.771995^2);

nE = sqrt(ne2);
end