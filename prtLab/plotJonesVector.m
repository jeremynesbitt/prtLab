% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function plotJonesVector(jin, w)
% Fancy plot of jones vector (with arrow!)
% Note for people like me who tend to forget.  
% The plots here are done by remembering that the E field varies
% in time as well as space, and these plot are over one period 
% of the time varying field, represented here by theta.
% V2 - if user asks for a ctr point + max size then scale accordingly
% this is added to provide 2d plots of vectors across some variable range
% w = [xo,yo, maxR]

assert(length(jin) == 2, 'Error!  Inproper jones matrix!');
%assert(isreal(jin(1)), 'Error!  Plots do not support imaginary first term')

magnitude = sqrt(sum(abs(jin).*abs(jin)));

if exist('w', 'var')
    magnitude = abs(w(3));
    jin = magnitude*jin;
end

theta = linspace(0,2*pi,101);
ex = exp(-1i*theta)*jin(1);
ey = exp(-1i*theta)*jin(2);

if exist('w','var')
    ex = ex + w(1);
    ey = ey + w(2);
end

%figJ = figure;
lineWidth = 6;
if exist('w','var')
    lineWidth = 2;
end
plot(real(ex),real(ey),'LineWidth',lineWidth, 'Color','r');

% Now the hard part.  drawing the arrow correctly.  
h = 0.1*magnitude; % height of arrow (relative to magnitude)
aT = 30*pi/180; % angle of arrow triangle 

if isreal(jin(1)) && isreal(jin(2)) % Linear pol
    phi = atan2(jin(2),jin(1));
    v1 = [jin(1);jin(2)]; 
    dr = magnitude*[cos(phi); sin(phi)];
    vb = v1-sqrt(2)*h*v1/norm(v1); % base should be along the vector v1
    dphi = magnitude*[-sin(phi);cos(phi)];
    v2 = vb+tan(aT)*dphi/4;
    v3 = vb-tan(aT)*dphi/4;
    triMat = [v1,v2,v3];
    if exist('w','var')
        triMat(1,:) = triMat(1,:) + w(1);
        triMat(2,:) = triMat(2,:) + w(2);
    end
    hold on;
    fill(triMat(1,:),triMat(2,:),'r');
    if ~exist('w','var')
        xlim([-magnitude,magnitude]);
        ylim([-magnitude,magnitude]);   
    end
else % some form of circular

    theta0 = -pi/6; % pick an angle.  TODO make optional input
    ex0 = exp(-1i*theta0)*jin(1);
    ey0 = exp(-1i*theta0)*jin(2);

    v1 = [real(ex0);real(ey0)];

    % Use finite differences to find dphi/dr.  May be a more elegant
    % way but this works for circular/elliptical
    dtta = .01;
    dex0 = exp(-1i*(theta0+dtta))*jin(1);
    dey0 = exp(-1i*(theta0+dtta))*jin(2); 
    dphi = 1./(theta0-dtta).*[real(ex0-dex0);real(ey0-dey0)];
    dphi = dphi*1/sqrt(sum(dphi.*dphi));
    dr = magnitude*[-dphi(2); dphi(1)];
        
    % Essentially switched vectors from linear
    %vb = v1-h*dphi;
    vb = v1-sqrt(2)*h*dphi;
    % Debug - plot dphi and dr.  
    %hold on;
    %plot([v1(1),vb(1)],[v1(2),vb(2)], 'r')
    %plot([v1(1),v1(1)-h*dr(1)],[v1(2),v1(2)-h*dr(2)],'g')        

        
    v2 = vb+tan(aT)*dr/4;
    v3 = vb-tan(aT)*dr/4;
    triMat = [v1,v2,v3];
    if exist('w','var')
        triMat(1,:) = triMat(1,:) + w(1);
        triMat(2,:) = triMat(2,:) + w(2);
    end    
    hold on;
    fill(triMat(1,:),triMat(2,:),'r');
    if ~exist('w','var')
        xlim([-magnitude,magnitude]);
        ylim([-magnitude,magnitude]);   
    end

end % Is real if

end % function
