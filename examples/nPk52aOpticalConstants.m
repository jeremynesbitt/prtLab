function n = nPk52aOpticalConstants(lambda)
%NPK52AOPTICALCONSTANTS Refractive index of Schott N-PK52A glass.
%
%   n = nPk52aOpticalConstants(lambda)
%
%   lambda is the vacuum wavelength in um. The coefficients are from the
%   Schott optical-glass overview workbook dated 13-Nov-2025. N-PK52A is
%   isotropic, so this function returns one refractive index.

B1 = 1.029607;
B2 = 0.1880506;
B3 = 0.736488165;
C1 = 0.00516800155;
C2 = 0.0166658798;
C3 = 138.964129;

lambda2 = lambda.^2;
n2 = 1 ...
    + (B1*lambda2)./(lambda2 - C1) ...
    + (B2*lambda2)./(lambda2 - C2) ...
    + (B3*lambda2)./(lambda2 - C3);

n = sqrt(n2);

end
