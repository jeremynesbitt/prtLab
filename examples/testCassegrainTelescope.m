addpath(genpath('..'));
T = exampleCassegrainTelescope();
k0 = [0; 0; 1];
x0 = [0; 4000;-6000];
x1 = [0; -4000;-6000];
d1 = polarizationRayTrace(T,k0,x0,[1;0]);
options = d1.options;
options.minAmplitude = .05;
d1 = polarizationRayTrace(T,k0,x0,[1;0]);
d2 = polarizationRayTrace(T,k0,x1,[1;0]);
Tplot = updateClearAperturesFromRayTrace(T, {d1,d2});
figure;plotPrtLensCrossSection(Tplot)
plotPrtRayTrace({d1,d2});

d = 4161.65; 
[X,Y] = meshgrid(linspace(-d,d,50)); % assume 10mm diameter pupil
[~,RHO] = cart2pol(X,Y);
ftr = RHO>=1000 & RHO<=d; % 
E_in = [1;0];
k_in = [0; 0; 1];
Jall = nan([size(X), 2, 2]);
OPL = nan(size(X));

coord = struct('type',"doublePole", ...
    'a_loc',k_in, ...
    'x_o',[1;0;0]);

for ii=1:size(X,1)
    for jj=1:size(X,2)
        pos_in = [X(ii,jj); Y(ii,jj); -6000];
        rayOutput = polarizationRayTrace(T, k_in, pos_in, E_in);
        if ~isempty(rayOutput.finalRayIds)
        finalRay = rayOutput.rays(rayOutput.finalRayIds(1));
        P_tot = finalRay.P;
        J = transformPtoJones(P_tot, k_in, finalRay.k, coord);
        Jall(ii,jj,:,:) = J;
        OPL(ii,jj) = finalRay.OPL;
        end
    end
end

figure;plotJonesPupil(Jall,ftr);
OPD = OPL - mean(OPL(ftr & isfinite(OPL)), 'all', 'omitnan');
clim = [min(OPD(ftr)), max(OPD(ftr))];
%figure;imagesc(OPD, 'AlphaData', ftr & isfinite(OPD)); axis image; colorbar;
figure;imagesc(OPD.*ftr, clim); axis image; colorbar;
title('OPD relative to pupil mean');

n_ideal = 1E15;
T.IndexData{2}.n = n_ideal;
T.IndexData{3}.n = n_ideal;

n_gold = 0.8 + 1.8*sqrt(-1);
nMirror = 0.958 + 6.69*1i;
T.IndexData{2}.n = n_gold;
T.IndexData{3}.n = n_gold;
