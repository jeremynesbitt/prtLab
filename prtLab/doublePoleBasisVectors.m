function [x_loc, y_loc] = doublePoleBasisVectors(k, a_loc, x_o)
% doublePoleBasisVectors - Compute double pole local coordinates on unit sphere
%
% [x_loc, y_loc] = doublePoleBasisVectors(k, a_loc, x_o) computes the double
% pole local coordinate basis vectors for propagation vectors k.
%
% Inputs:
%   k     - 3x1 array of unit propagation vectors
%   a_loc - 3x1 axis vector (anti-pole direction, singularity at -a_loc)
%   x_o   - 3x1 reference x-basis vector at the anti-pole (must be perp to a_loc)
%
% Outputs:
%   x_loc - Nx3 array of local x-basis vectors
%   y_loc - Nx3 array of local y-basis vectors
%
% Based on Eq 11.10 of Chipman, Lam, Young:
%   (x_Loc, y_Loc) = (R * x_o, R * (a_Loc x x_o))
% where R is a rotation about r = k x a_Loc by angle theta = -acos(k . a_Loc).
%
% Uses Rodrigues' rotation formula. Singular at k = -a_loc (the double pole).

% Original design assumed 1x3 but all other code assumes 3x1 so just fixing
% with a transpose here;
k=k';
a_loc=a_loc';
x_o = x_o';

a = a_loc(:)' / norm(a_loc);  % ensure row, normalized
x_o = x_o(:)' / norm(x_o);   % ensure row, normalized

% Verify x_o is perpendicular to a_loc
assert(abs(dot(a, x_o)) < 1e-6, 'x_o must be perpendicular to a_loc');

y_o = cross(a, x_o);  % reference y at anti-pole

N = size(k, 1);
x_loc = zeros(N, 3);
y_loc = zeros(N, 3);

for ii = 1:N
    ki = k(ii, :);
    cos_theta = dot(ki, a);

    if cos_theta > 1 - 1e-10
        % k is at the anti-pole (k = a_loc): no rotation needed
        x_loc(ii, :) = x_o;
        y_loc(ii, :) = y_o;
    elseif cos_theta < -1 + 1e-10
        % k is at the double pole (k = -a_loc): singular
        % Set to arbitrary basis (similar to dipole pole handling)
        x_loc(ii, :) = x_o;
        y_loc(ii, :) = y_o;
    else
        % General case: rotate x_o and y_o using Rodrigues' formula
        % Rotation axis: r = k x a (unnormalized)
        r = cross(ki, a);
        r = r / norm(r);  % unit rotation axis

        % Rotation angle: theta = -acos(k . a)
        theta = -acos(cos_theta);
        ct = cos(theta);
        st = sin(theta);

        % Rodrigues: v_rot = v*cos(t) + (r x v)*sin(t) + r*(r.v)*(1-cos(t))
        x_loc(ii, :) = x_o * ct + cross(r, x_o) * st + r * dot(r, x_o) * (1 - ct);
        y_loc(ii, :) = y_o * ct + cross(r, y_o) * st + r * dot(r, y_o) * (1 - ct);
    end
end

end
