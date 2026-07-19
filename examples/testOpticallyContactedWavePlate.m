addpath(genpath('..'));

T = exampleAchromaticWavePlateSystem();
lambda = T.Properties.UserData.lambda ;
k = [0; sin(20*pi/180); cos(20*pi/180)];
rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1]);
options = rayOutput.options;
options.minAmplitude = .01;
rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

OPL_in_waves = OPL/lambda

% Look at wavelength dependence
lambda = 0.633; % um
no2 = 1 ...
    + 0.28604141 ...
    + (1.07044083*lambda.^2)./(lambda.^2 - 1.00585997e-2) ...
    + (1.10202242*lambda.^2)./(lambda.^2 - 100);

nO = sqrt(no2);

ne2 = 1 ...
    + 0.28851804 ...
    + (1.09509924*lambda.^2)./(lambda.^2 - 1.02101864e-2) ...
    + (1.15662475*lambda.^2)./(lambda.^2 - 100);

nE = sqrt(ne2);

tQwp = 1.25 * lambda / (nO - nE) / 2;
opticAxis_1 = [-sqrt(2)/2; -sqrt(2)/2; 0];
opticAxis_2 = [sqrt(2)/2; sqrt(2)/2; 0];

T = modAchromaticWavePlateSystem(lambda,nO,nE,nO,nE,tQwp,tQwp,opticAxis_1,opticAxis_2);
rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

OPL_in_waves = OPL/lambda

% Test vs wavelength
wlArray = 1/1000*linspace(400,700,101);

for ii=1:length(wlArray)
    lambda = wlArray(ii);
    % Calcite vs wavelength
    no2 = 1 ...
        + 0.28604141 ...
        + (1.07044083*lambda.^2)./(lambda.^2 - 1.00585997e-2) ...
        + (1.10202242*lambda.^2)./(lambda.^2 - 100);
    
    nO = sqrt(no2);
    
    ne2 = 1 ...
        + 0.28851804 ...
        + (1.09509924*lambda.^2)./(lambda.^2 - 1.02101864e-2) ...
        + (1.15662475*lambda.^2)./(lambda.^2 - 100);
    
    nE = sqrt(ne2);

    T = modAchromaticWavePlateSystem(lambda,nO,nE,nO,nE,tQwp,tQwp,opticAxis_1,opticAxis_2);
    rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

    OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

    OPL_AchromaticWaveplate(ii,1) = OPL/lambda;
end


%% Test output
E_in = [1;0];
[X,Y] = meshgrid(linspace(-5,5,17)); % assume 10mm diameter pupil
[tX, tY] = meshgrid(linspace(-5,5,17)); % 

coord = struct('type',"doublePole", ...
    'a_loc',[0;0;1], ...
    'x_o',[1;0;0]);

for ii=1:size(tX,1)
    for jj=1:size(tX,2)
        k_in = calcKFromThetaXThetaY(tX(ii,jj), tY(ii,jj));
        pos_in = [X(ii,jj); Y(ii,jj); 0];
        rayOutput = polarizationRayTrace(T, k_in, pos_in, E_in);
        finalRays = rayOutput.finalRayIds;
        P_tot = zeros(3,3);
        for kk=1:length(finalRays)
            P_tot = P_tot+rayOutput.rays(finalRays(kk)).P;
        end

        J = transformPtoJones(P_tot, k_in, rayOutput.rays(rayOutput.finalRayIds(1)).k, coord);
        Jall(ii,jj,:,:) = J;
      
    end
end

plotPolarizationEllipsesAcrossField(X,Y,Jall, E_in); title('True 0 Order Plate');
plotPolarizationEllipsesAcrossField(X,Y,Jall_m, E_in); title('Multi Order plate');


% Calcite vs wavelength
no2 = 1 ...
    + 0.28604141 ...
    + (1.07044083*lambda.^2)./(lambda.^2 - 1.00585997e-2) ...
    + (1.10202242*lambda.^2)./(lambda.^2 - 100);

no = sqrt(no2);

ne2 = 1 ...
    + 0.28851804 ...
    + (1.09509924*lambda.^2)./(lambda.^2 - 1.02101864e-2) ...
    + (1.15662475*lambda.^2)./(lambda.^2 - 100);

ne = sqrt(ne2);