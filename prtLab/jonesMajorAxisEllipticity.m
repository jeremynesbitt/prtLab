function [psi, ellipticity, majorAxis, a, b] = jonesMajorAxisEllipticity(E)
% JONESMAJORAXISELLIPTICITY  Major axis orientation and ellipticity.
%
%   [psi, ellipticity] = jonesMajorAxisEllipticity(E) computes the major
%   axis orientation psi and ellipticity epsilon = b/a from a Jones vector
%   E = [Ex; Ey] = [Ax*exp(-1i*phix); Ay*exp(-1i*phiy)].
%
%   [psi, ellipticity, majorAxis, a, b] also returns the unit vector along
%   the major axis, the semi-major axis length a, and the semi-minor axis
%   length b.
%
%   E may be a 2xN matrix; outputs are 1xN except majorAxis, which is 2xN.
%
%   Definitions used:
%     phi = phiy - phix
%     tan(2*psi) = 2*Ax*Ay*cos(phi) / (Ax^2 - Ay^2)
%     a = sqrt(Ax^2*cos(psi)^2 + Ay^2*sin(psi)^2 ...
%              + 2*Ax*Ay*cos(psi)*sin(psi)*cos(phi))
%     b = sqrt(Ax^2*sin(psi)^2 + Ay^2*cos(psi)^2 ...
%              - 2*Ax*Ay*cos(psi)*sin(psi)*cos(phi))
%     epsilon = b/a
%
%   See also: jonesToStokes, stokesToJones

if size(E, 1) ~= 2
    error('jonesMajorAxisEllipticity: E must be a 2xN matrix (2 rows).');
end

Ex = E(1, :);
Ey = E(2, :);

Ax = abs(Ex);
Ay = abs(Ey);

% With E = A*exp(-1i*phase), phi = phiy - phix = angle(Ex) - angle(Ey).
cosPhi = ones(size(Ax));
nonzero = Ax > 0 & Ay > 0;
cosPhi(nonzero) = real(Ex(nonzero) .* conj(Ey(nonzero))) ./ ...
                  (Ax(nonzero) .* Ay(nonzero));
cosPhi = min(1, max(-1, cosPhi));

psi = 0.5 * atan2(2 * Ax .* Ay .* cosPhi, Ax.^2 - Ay.^2);

cosPsi = cos(psi);
sinPsi = sin(psi);

aSquared = Ax.^2 .* cosPsi.^2 + Ay.^2 .* sinPsi.^2 + ...
           2 * Ax .* Ay .* cosPsi .* sinPsi .* cosPhi;
bSquared = Ax.^2 .* sinPsi.^2 + Ay.^2 .* cosPsi.^2 - ...
           2 * Ax .* Ay .* cosPsi .* sinPsi .* cosPhi;

a = sqrt(max(0, aSquared));
b = sqrt(max(0, bSquared));

ellipticity = zeros(size(a));
nonzeroMajor = a > 0;
ellipticity(nonzeroMajor) = b(nonzeroMajor) ./ a(nonzeroMajor);

majorAxis = [cosPsi; sinPsi];
end
