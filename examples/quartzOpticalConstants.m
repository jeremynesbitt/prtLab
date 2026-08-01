% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [nO, nE] = quartzOpticalConstants(lambda)
%lambda is in um
no2 = 1 ...
    + 0.28604141 ...
    + (1.07044083*lambda.^2)./(lambda.^2 - 1.00585997e-2) ...
    + (1.10202242*lambda.^2)./(lambda.^2 - 100);

nO = sqrt(no2);

ne2 = 1 ...
    + 0.28851804 ...
    + (1.09509924*lambda.^2)./(lambda.^2 - 1.02101864e-2) ...
    + (1.15662475*lambda.^2)./(lambda.^2 - 100);

nE = sqrt(ne2);
end