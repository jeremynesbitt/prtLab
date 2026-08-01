% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function plotPolarizationEllipsesAcrossField(X, Y, Jout, Ein, ftr)
%PLOTPOLARIZATIONELLIPSESACROSSFIELD Plot output polarization ellipses.
%
%   Jout may be either:
%       M x N x 2 x 2   preferred, stores each Jones matrix directly
%       M x N x 4       legacy order [J11 J12 J21 J22]

figure;
mag = (Y(2)-Y(1))/2;

if ~exist('ftr','var')
    [~, RHO] = cart2pol(X,Y);
    ftr = RHO <= sqrt(max(X(:)).^2 + max(Y(:)).^2);
end
for ii=1:size(X,1)
    for jj=1:size(X,2)
        if ftr(ii,jj) == 1
            J = localJonesMatrix(Jout, ii, jj);
            tmp = J * Ein;
            plotJonesVector(tmp, [-X(ii,jj), Y(ii,jj), mag]);
        end
    end
end
title('Output Polarization');

end

function J = localJonesMatrix(Jout, ii, jj)
if ndims(Jout) == 4 && size(Jout,3) == 2 && size(Jout,4) == 2
    J = squeeze(Jout(ii,jj,:,:));
elseif ndims(Jout) == 3 && size(Jout,3) == 4
    J = [Jout(ii,jj,1), Jout(ii,jj,2); ...
         Jout(ii,jj,3), Jout(ii,jj,4)];
else
    error('plotPolarizationEllipsesAcrossField:InvalidJonesArray', ...
        'Jout must be MxNx2x2 or MxNx4.');
end
end
