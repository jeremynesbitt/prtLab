
addpath(genpath('..'));

lambda = 0.589; % um
nGlass = 1.4965; %1.495;

rhomb = createFresnelRhombSolid(nGlass, ...
    Length=.5*28, ...
    Height=8, ...
    Width=8, ...
    Shear=.5*22.0417, ...
    lambda=lambda, ...
    lambdaUnits="um");

coord = struct('type',"doublePole", ...
    'a_loc',[0;0;1], ...
    'x_o',[1;0;0]);

options = prtDefaultOptions();
options.maxInteractions = 8;
options.minAmplitude = 1e-4;

tta = 0;
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

disp("Hit sequence:");
for ii = 1:numel(rayOutput.interactions)
    data = rayOutput.interactions(ii).interceptData;
    fprintf("  %s / %s (%s)\n", ...
        data.solidName, data.faceName, data.direction);
end

figure;
plotPrtSystem3D(rhomb, rayOutput, ...
    RaySelection="dominant", ...
    PostExtend=4);
title("Fresnel rhomb polarization evolution");


angles = 0:5;

for ii=1:length(angles)
tta = angles(ii);
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

    finalRay = selectDominantFinalRay(rayOutput);
    P = finalRay.P;
    Q = finalRay.Q;
    kOut = finalRay.k;
    Jm = transformPtoJones(P, kIn, kOut, coord);


    % Compute Jout via rotation matrix
    ttaRot = acos(dot(kIn,[0;0;-1]));
    kx = kIn(1);
    ky = kIn(2);
    kz = kIn(3);
    H = kx^2+ky^2;

    if H==0
        U = eye(3);
    else
        U = 1/H*[kx^2*cos(ttaRot)+ky^2, kx*ky*(cos(ttaRot)-1), -sqrt(H)*kx*sin(ttaRot) ;...
              kx*ky*(cos(ttaRot)-1), kx^2+ky^2*cos(ttaRot), -sqrt(H)*ky*sin(ttaRot) ;...
              sqrt(H)*kx*sin(ttaRot), sqrt(H)*ky*sin(ttaRot), H*cos(ttaRot)];    
    end
    Eout = inv(U)*finalRay.fieldE;

retE(ii) = angle(Eout(1))-angle(Eout(2));
retJ(ii) = calcRetardanceFromJonesMatrix(Jm);
retQP(ii) = calcRetardanceByQPMethod(Q,P);
end

% Test vs wavelength
wls = 1/1000*linspace(400,700,21);
tta = 0;
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

for ii=1:length(wls)

lambda = wls(ii); % um
nGlass = nPk52aOpticalConstants(lambda);

rhomb = createFresnelRhombSolid(nGlass, ...
    Length=.5*28, ...
    Height=8, ...
    Width=8, ...
    Shear=.5*22.0417, ...
    lambda=lambda, ...
    lambdaUnits="um");

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

finalRay = selectDominantFinalRay(rayOutput);
P = finalRay.P;
Q = finalRay.Q;
kOut = finalRay.k;
Jm = transformPtoJones(P, kIn, kOut, coord);

retardanceVsWavelength(ii) = calcRetardanceFromJonesMatrix(Jm);

end

% Okay let's do a design study.  Step 1: find the rhomb geometry s.t. the
% retardance at the first interface is pi/4.

% Test vs wavelength
shears = linspace(10,12,101);
tta = 0;
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

lambda = .589; % um
nGlass = nPk52aOpticalConstants(lambda);

for ii=1:length(shears)

rhomb = createFresnelRhombSolid(nGlass, ...
    Length=14, ...
    Height=8, ...
    Width=8, ...
    Shear=shears(ii), ...
    lambda=lambda, ...
    lambdaUnits="um");

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

finalRay = selectDominantFinalRay(rayOutput);
P = finalRay.P;
kOut = finalRay.k;
Jm = transformPtoJones(P, kIn, kOut, coord);

retardanceVsShear(ii) = calcRetardanceFromJonesMatrix(Jm);

end


shear = 10.84; % eyeball fit

% Now do vs wavelength

% Test vs wavelength
wls = 1/1000*linspace(400,700,21);
tta = 0;
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

for ii=1:length(wls)

lambda = wls(ii); % um
nGlass = nPk52aOpticalConstants(lambda);

rhomb = createFresnelRhombSolid(nGlass, ...
    Length=14, ...
    Height=8, ...
    Width=8, ...
    Shear=shear, ...
    lambda=lambda, ...
    lambdaUnits="um");

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

finalRay = selectDominantFinalRay(rayOutput);
P = finalRay.P;
Q = finalRay.Q;
kOut = finalRay.k;
Jm = transformPtoJones(P, kIn, kOut, coord);

retardanceVsWavelength(ii) = calcRetardanceFromJonesMatrix(Jm);

end
