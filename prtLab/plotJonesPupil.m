% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function plotJonesPupil(Jmat, ftr)
% Create two figures, one for real part and one for imaginary part
% Jmat should be nxnx2x2

if ~exist('ftr', 'var')
    % default scale set within +/- 1
    [X,Y] = meshgrid(linspace(-1,1,size(Jmat,1)));
    [~, RHO] = cart2pol(X,Y);
    ftr = RHO<1;
end

hR = figure;
tR = tiledlayout(2,2,'TileSpacing','Compact');
hI = figure;
tI = tiledlayout(2,2,'TileSpacing','Compact');

for ii=1:2
    for jj=1:2
        J = Jmat(:,:,ii,jj);
        
        % For the phase, to avoid phase jumps when the sign of the real
        % part changes use the convention below.  Note that this does not
        % change the meaning of J. 

        % If J is small then the phase will not be that meaningful
        mask = ftr & abs(J) > 1e-4*max(abs(J(ftr)));


        phi0 = angle(mean(J(mask)));
        %phi0 = 0; % can choose a different mean value
        Jrot = J * exp(-1i*phi0);

        % Keep track of the sign of J
        sgn = sign(real(Jrot));
        sgn(sgn == 0) = 1;
        

        A_signed = real(Jrot);
        % This will guarantee the phase does not flip sign as the real part
        % changes sign.
        Jtmp = angle(Jrot .* sgn); 
        % For plotting purposes, nan the values that are too small.
        Jtmp(~mask) = nan;

        nexttile(tR); 
        clim = [min(A_signed(ftr)), max(A_signed(ftr))];
        imagesc(flipud(A_signed).*ftr, clim); colorbar;

        nexttile(tI); 
        if (sum(~isnan(Jtmp(:)))) > 0
        clim = [min(Jtmp(ftr)), max(Jtmp(ftr))];
        if clim(1)~=clim(2)
        imagesc(flipud(Jtmp).*ftr, clim); colorbar;
        else
        imagesc(flipud(Jtmp).*ftr); colorbar;    
        end
    end
end

% 3. Add an overall title
figTitle = title(tR, 'Real Part of Jones Pupil');

% (Optional) Customize the title
figTitle.FontSize = 16;
figTitle.FontWeight = 'bold';

% 3. Add an overall title
figTitle = title(tI, 'Imaginary Part of Jones Pupil [rad]');

% (Optional) Customize the title
figTitle.FontSize = 16;
figTitle.FontWeight = 'bold';

end % function