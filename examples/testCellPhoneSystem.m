% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

addpath(genpath('..'));

%T = exampleCellPhoneLensEx9System();
T = exampleCellPhoneLensEx3System();

k = [0; sin(20*pi/180); cos(20*pi/180)];
x1 = [0; 0.776; -0.16];
x2 = [0; -0.776; -0.16];

d1 = polarizationRayTrace(T, k, x1, [0;1]);
d2 = polarizationRayTrace(T, k, x2, [0;1]);

Tplot = updateClearAperturesFromRayTrace(T, {d1, d2}, ...
    'Margin', 1.15, ...
    'Minimum', 0.25);

T = Tplot;

figure;
plotPrtLensCrossSection(T);
plotPrtRayTrace({d1, d2});


rayOutput = polarizationRayTrace(T, [.0858;.173;.9117], [0;0;0], [0;1]);
rayOutput2 = polarizationRayTrace(T, [-.086;-.173;.981], [-.086;-.173;-.981], [0;1]);
%figure;
%plotPrtLensCrossSection(T);
%plotPrtRayTrace(rayOutput);

figure;
plotPrtLensCrossSection(T);
plotPrtRayTrace([rayOutput;rayOutput2]);

% Test on axis 

x = [0;  0.741676704231623; -0.1];
k = [0;  0;                 1];

rayOutput = polarizationRayTrace(T, k, x, [0;1]);
x = [0;  -0.741676704231623; -0.16];
k = [0;  0;                 1];
rayOutput2 = polarizationRayTrace(T, k, x, [0;1]);
figure;
plotPrtLensCrossSection(T);
plotPrtRayTrace([rayOutput;rayOutput2]);

%k = [0; sin(15*pi/180); cos(15*pi/180)];
k = [0; sin(7*pi/180); cos(7*pi/180)];
%k = [0; sin(20*pi/180); cos(20*pi/180)];
%k = [ 0.0000000000	 ;  0.3451737095	;   0.9385388166];

%[X,Y] = meshgrid(linspace(-0.995,0.995,17));

[X,Y] = meshgrid(linspace(-0.741676704231623,0.741676704231623,17));

Ein = [1;0]; 
%[~, sIn] = calcPandSUnitVectors(k, )

%% Look at rays for debug
ii=9;jj=9; rayOutput = polarizationRayTrace(T, k, [X(ii,jj) ; Y(ii,jj) ; -0.16], Ein);
ii=17;jj=9; rayOutput2 = polarizationRayTrace(T, k, [X(ii,jj) ; Y(ii,jj) ; 0], Ein);
figure;
plotPrtLensCrossSection(T);
plotPrtRayTrace([rayOutput;rayOutput2]);

x_loc = [];
y_loc = [];
k_loc = [];

[~, RHO] = cart2pol(X,Y);
ftr = RHO <= 0.741676704231623;

for ii=1:size(X,1)
    for jj=1:size(X,2)
        rayOutput = polarizationRayTrace(T, k, [X(ii,jj) ; Y(ii,jj) ; 0], Ein);
        if isempty(rayOutput.finalRayIds)
            Pall(ii,jj,1:9) = 0;
        else
        finalRay = rayOutput.rays(rayOutput.finalRayIds(1));
        Ptotal = finalRay.P;
        Pall(ii,jj,:) = Ptotal(:);

        % Look at output in double pole coordinates.
       
         % TODO:  Define Ein as s pol.  Compute Pout = Ptot*Ein, then 
         % double pol xform Pout for plotting
        x_o = [1;0;0];
        [x_dp, y_dp] = doublePoleBasisVectors(k, k, x_o);
        U3 = [x_dp', y_dp', k];
        [x_dp, y_dp] = doublePoleBasisVectors(finalRay.k, k, x_o);
        U2 = [x_dp', y_dp', finalRay.k];
        Jtmp = U2\(Ptotal*U3);
        Jout(ii,jj,1:4) = [Jtmp(1,1), Jtmp(1,2), Jtmp(2,1), Jtmp(2,2)];



        % Look at output in s/p coordinates
        [p_i, s_i] = calcPandSUnitVectors(k,rayOutput.interactions(1).normal);
        U3 = [s_i, p_i, k];

        idx = 11;  % or whatever surface you mean
        k_o = rayOutput.interactions(idx).children(1).k;     % transmitted direction after surface idx
        n_o = rayOutput.interactions(idx).normal;
        [p_o, s_o] = calcPandSUnitVectors(k_o, n_o);
        %[p_o, s_o] = calcPandSUnitVectors(finalRay.k, rayOutput.interactions(end).normal);
        Ptotal = rayOutput.rays(idx+1).P;
        Pall(ii,jj,:) = Ptotal(:);
        U2 = [s_o, p_o, k_o];
         Jtmp = U2\(Ptotal*U3);
        Jout(ii,jj,1:4) = [Jtmp(1,1), Jtmp(1,2), Jtmp(2,1), Jtmp(2,2)];

        % Test for singularity
        if isreal(finalRay.k) && ftr(ii,jj) == 1
        x_loc = [x_loc ; s_o'];
        y_loc = [y_loc ; p_o'];
        k_loc = [k_loc ; k_o'];
        end
            
        end
        %
        %Pall(ii,jj,:) = tmpOut.interactions(end).P.transmitted(:);

    end
end



figure;plotPMatrices(Pall.*ftr);

figure;plotJonesPupil(Jout.*ftr)


figure;
mag = (Y(2)-Y(1))/2;
for ii=1:size(X,1)
    for jj=1:size(X,2)
        if ftr(ii,jj) == 1
        plotJonesVector(Ein, [-X(ii,jj), Y(ii,jj), mag]);
        end
    end
end
title('Input Polarization');

figure;
mag = (Y(2)-Y(1))/2;
for ii=1:size(X,1)
    for jj=1:size(X,2)
        if ftr(ii,jj) == 1
        tmp = squeeze([Jout(ii,jj,1), Jout(ii,jj,2) ; Jout(ii,jj,3), Jout(ii,jj,4)] *Ein);
        plotJonesVector(tmp, [-X(ii,jj), Y(ii,jj), mag]);
        end
    end
end
title('Output Polarization');

ii=9;jj=9; rayOutput = polarizationRayTrace(T, k, [X(ii,jj) ; Y(ii,jj) ; -0.16], Ein);
k_chief = rayOutput.rays(end).k ;
tY = asin(k_chief(2))*180/pi;
opts_k = struct('showSphere', false, 'arrowScale', 0.02, ...
    'sphereAlpha', 0.1, 'lineWidth', 1, 'viewAngle', [0, tY, -90]);
figure;plotBasisVectorsOnSphere(k_loc,x_loc, [], opts_k);
