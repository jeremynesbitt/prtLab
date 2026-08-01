% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function [v_rot, Rz] = rotateToZeroX(v)
%ROTATETOZEROX Rotate a 3D vector so that its x-component becomes zero.
%
%   v_rot = ROTATETOZEROX(v)
%
%   Input:
%       v     - 3×1 or 1×3 vector [x y z]
%
%   Output:
%       v_rot - rotated vector with v_rot(1) = 0

    % Ensure column vector
    v = v(:);

    x = v(1);
    y = v(2);
    z = v(3);

    % If x is already zero, nothing to do
    if abs(x) < 1e-14
        v_rot = v;
        Rz = eye(3);
        return
    end

    % Rotation angle around z so that x' = 0
    % tan(phi) = -x / y
   phi = atan2(x, y);   % gives correct quadrant and avoids division by zero

    % Rotation matrix about z-axis
    Rz = [ cos(phi)  -sin(phi)   0;
           sin(phi)   cos(phi)   0;
                0          0     1 ];

    % Rotate vector
    v_rot = Rz * v ;

    % if sign(v_rot(2)) ~= sign(x)
    %     Rz = inv(Rz);
    %     v_rot = Rz*v;
    % end

end
