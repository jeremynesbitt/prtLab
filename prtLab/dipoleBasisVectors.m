function [x_loc, y_loc] = dipoleBasisVectors(k, a_loc)
% dipoleBasisVectors - Compute dipole local coordinates on unit sphere
%
% [x_loc, y_loc] = dipoleBasisVectors(k, a_loc) computes the dipole
% (latitude/longitude) local coordinate basis vectors for propagation
% vectors k with axis vector a_loc.
%
% Inputs:
%   k     - Nx3 array of unit propagation vectors
%   a_loc - 1x3 axis vector (defines pole locations at +/- a_loc)
%
% Outputs:
%   x_loc - Nx3 array of local x-basis vectors (along lines of constant latitude)
%   y_loc - Nx3 array of local y-basis vectors (along lines of constant longitude)
%
% Based on Eq 11.5 of Chipman, Lam, Young:
%   x_loc = (a x k) / |a x k|
%   y_loc = (k x a x k) / |a x k|
%
% Singular at k = +/- a_loc (the poles). At these points, basis vectors
% are set to a default orthonormal pair.

a = a_loc(:)' / norm(a_loc);  % ensure row vector, normalized
N = size(k, 1);

x_loc = zeros(N, 3);
y_loc = zeros(N, 3);

% Compute a x k for all points
axk = cross(repmat(a, N, 1), k, 2);  % Nx3
mag = sqrt(sum(axk.^2, 2));          % Nx1

% Threshold for singularity (k near +/- a_loc)
tol = 1e-6;
good = mag > tol;

% Regular points: Eq 11.5
x_loc(good, :) = axk(good, :) ./ mag(good);
% y_loc = k x x_loc
y_loc(good, :) = cross(k(good, :), x_loc(good, :), 2);

% Singular points: pick arbitrary orthonormal basis in transverse plane
if any(~good)
    idx = find(~good);
    for ii = 1:length(idx)
        ki = k(idx(ii), :);
        % Pick a reference vector not parallel to ki
        if abs(ki(1)) < 0.9
            ref = [1, 0, 0];
        else
            ref = [0, 1, 0];
        end
        b1 = cross(ki, ref);
        b1 = b1 / norm(b1);
        b2 = cross(ki, b1);
        x_loc(idx(ii), :) = b1;
        y_loc(idx(ii), :) = b2;
    end
end

end
