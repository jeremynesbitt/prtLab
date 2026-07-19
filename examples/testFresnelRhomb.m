
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

angles = 0:5;

for ii=1:length(angles)
tta = angles(ii);
kIn = [0; sin(tta*pi/180); cos(tta*pi/180)];
xIn = [0; 5; -2];
Ein = [1; 1]/sqrt(2);

rayOutput = polarizationRayTraceSolid(rhomb, kIn, xIn, Ein, options);

P = rayOutput.rays(rayOutput.finalRayIds(1)).P;
Q = rayOutput.rays(rayOutput.finalRayIds(1)).Q;
kOut = rayOutput.rays(rayOutput.finalRayIds(1)).k;
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
    Eout = inv(U)*rayOutput.rays(rayOutput.finalRayIds(1)).fieldE;

retE(ii) = angle(Eout(1))-angle(Eout(2));
retJ(ii) = calcRetardanceFromJonesMatrix(Jm);
retQP(ii) = calcRetardanceByQPMethod(Q,P);
end

