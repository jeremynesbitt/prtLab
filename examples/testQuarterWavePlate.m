addpath(genpath('..'));

T = exampleQuarterWavePlateSystem();
lambda = T.Properties.UserData.lambda ;
k = [0; 0; 1];
rayOutput = polarizationRayTrace(T, k, [0;0;0], [0;1]);

OPL = rayOutput.rays(rayOutput.finalRayIds(1)).OPL - rayOutput.rays(rayOutput.finalRayIds(2)).OPL; 

OPL_in_waves = OPL/lambda

finalRays = rayOutput.finalRayIds;
P_tot = zeros(3,3);
for kk=1:length(finalRays)
    P_tot = P_tot+rayOutput.rays(finalRays(kk)).P;
end
Eout_x = P_tot*[1;0;0];
Eout_y = P_tot*[0;1;0]