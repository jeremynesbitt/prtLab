% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function retardance = calcRetardanceFromJonesMatrix(J)
% From PLAOS Chapter 17, equation 17.31

%   retardance = 2*acos( abs( tr(J) + det(J)/abs(det(J))*tr(J') ) ...
%                   / ( 2*sqrt(tr(J'*J) + 2*abs(det(J))) ) );

    d = det(J);

    retardance = 2*acos( ...
        abs(trace(J) + (d/abs(d))*trace(J')) / ...
        (2*sqrt(trace(J'*J) + 2*abs(d))) ...
        );
end

