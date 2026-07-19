addpath(genpath('..'));

T = exampleAchromaticWavePlateSystem();
t1 = T.Thickness(2);
tgap = T.Thickness(3);
t2 = T.Thickness(4);
lambda = T.Properties.UserData.lambda ;
opticAxis_1 = T.AxisData{2}.opticAxis;
opticAxis_2 = T.AxisData{4}.opticAxis;

% Axes should be crossed
assert(dot(opticAxis_1,opticAxis_2)==0, 'Optical Axes Should be crossed');

%opticAxis_1 = [1;0;0];
%opticAxis_2 = [0;1;0];
k = [0; 0; 1];
rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1]);
options = rayOutput.options;
options.minAmplitude = .01;
rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

OPL_in_waves = OPL/lambda

finalRays = rayOutput.finalRayIds;
P_tot = zeros(3,3);
for kk=1:length(finalRays)
    P_tot = P_tot+rayOutput.rays(finalRays(kk)).P;
end

J = P_tot(1:2,1:2); % normal incidence okay

ev = eig(J)
ret = angle(ev(2)/ev(1)) / (2*pi);
ret = mod(ret, 1);

%k = [0; sin(20*pi/180); cos(20*pi/180)];

%T = modAchromaticWavePlateSystem(lambda,nO,nE,nO,nE,tQwp,tQwp,opticAxis_1,opticAxis_2);
%rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

%OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

%OPL_in_waves = OPL/lambda

% Test vs wavelength
wlArray = 1/1000*linspace(400,700,101);

for ii=1:length(wlArray)
    lambda = wlArray(ii);
    [nO1, nE1] = quartzOpticalConstants(lambda);
    [nO2, nE2] = mgF2OpticalConstants(lambda);

    T = modAchromaticWavePlateSystem(lambda,nO1,nE1,nO2,nE2,t1,tgap,t2,opticAxis_1,opticAxis_2);
    rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

        finalRays = rayOutput.finalRayIds;
        P_tot = zeros(3,3);
        for kk=1:length(finalRays)
            P_tot = P_tot+rayOutput.rays(finalRays(kk)).P;
        end
    Eout_x = P_tot*[1;0;0];
    Eout_y = P_tot*[0;1;0];

    % Second way - via J
    J = P_tot(1:2,1:2); % normal incidence okay
    
    ev = eig(J);
    ret = angle(ev(2)/ev(1)) / (2*pi);
    ret = mod(ret, 1); 

    ret = min(ret, 1-ret);
    retInWaves(ii,1) = ret;

    [psi, ellipticity, majorAxis, a, b] = jonesMajorAxisEllipticity(Eout_x(1:2));
    ellipseParams(ii,:) = [psi, ellipticity, a, b];
    OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

    OPL_AchromaticWaveplate(ii,1) = OPL/lambda;
    OPL_firstOrder(ii,1) = 1/lambda*((nO1-nE1)*t1 + (nO2-nE2)*t2);
    OPL_firstOrder_mm(ii,1) = ((nO1-nE1)*t1 + (nO2-nE2)*t2);
    opticalConstantsvsLambda(ii,:) = [nO1, nE1, nO2, nE2];
end

% Look at polarization glyph at different wavelengths
figure;
wlArray = [.400, .550, .700];
for ii=1:length(wlArray)
    lambda = wlArray(ii);
    [nO1, nE1] = quartzOpticalConstants(lambda);
    [nO2, nE2] = mgF2OpticalConstants(lambda);

    T = modAchromaticWavePlateSystem(lambda,nO1,nE1,nO2,nE2,t1,tgap,t2,opticAxis_1,opticAxis_2);
    rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1], options);

        finalRays = rayOutput.finalRayIds;
        P_tot = zeros(3,3);
        for kk=1:length(finalRays)
            P_tot = P_tot+rayOutput.rays(finalRays(kk)).P;
        end
    Eout = P_tot*[1;0;0] ;      
    plotJonesVector(Eout(1:2)); % okay for normal incidence
end

