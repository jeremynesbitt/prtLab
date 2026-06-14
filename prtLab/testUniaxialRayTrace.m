T = exampleQuarterWavePlateSystem();

rayOutput = polarizationRayTrace(T, [0;0;1], [0;0;0], [0;1]);

% Compare with old way
s.t = T.Thickness(2);
s.lambda = T.Properties.UserData.lambda;
s.thetaX = 0;
s.thetaY = 0;
s.ada = [0;0;1];
opticAxis = T.AxisData{2}.opticAxis;
s.cosX = opticAxis(1);
s.cosY = opticAxis(2);
s.cosZ = opticAxis(3);
s.nO = T.IndexData{2}.nO;
s.nE = T.IndexData{2}.nE;

sout = birefringentRayTracePlanarInterface(s);
soutP = birefringentToIsotropicPMatrix(sout);
