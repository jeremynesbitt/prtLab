function kinc = calcKFromThetaXThetaY(thetaX,thetaY)
% Assume thetaX, thetaY are in degrees
deg = pi/180;

% Following convention in chapter 11 of PLAOS
kincNorm = abs(sqrt(tan(thetaX*deg)^2+tan(thetaY*deg)^2+1));
kinc = 1/kincNorm*[tan(thetaX*deg); tan(thetaY*deg); 1];


end