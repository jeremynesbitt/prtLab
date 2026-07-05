function interaction = traceUniaxialToUniaxial(interaction, mediumIn, mediumOut, ray, hit, normal, options)
%TRACEUNIAXIALTOUNIAXIAL trace rays from uniaxial to uniaxial materials
%Based on PLAOS Section 19.7.4


deg = pi/180;

k_inc = interaction.incident.k;
S_inc = interaction.incident.S;
E_inc = interaction.incident.modeE;
H_inc = interaction.incident.modeH;
iStruct = mediumIn.IndexData{1};
nO_i = iStruct.nO;
nE_i = iStruct.nE;
oStruct = mediumOut.IndexData{1};
nO_t = oStruct.nO;
nE_t = oStruct.nE;
if strcmp(interaction.incident.mode, "ordinary")
    n_inc = iStruct.nO;
else
    n_inc = ray.metadata.n;
end
ada = interaction.normal;

axisData = mediumOut.AxisData{1};
opticAxis_t = axisData.opticAxis;
if norm(opticAxis_t) ~= 1
    opticAxis_t = opticAxis_t/norm(opticAxis_t);
end


Iin  = dot(cross(E_inc, conj(H_inc)), ada);
IinField = dot(cross(ray.fieldE, conj(ray.fieldH)), ada);

E_n = cross(S_inc,E_inc);
E_n = E_n/norm(E_n);



% The analytic solutions below require alpha to be in the yz plane.
% rotate the coordinate system if needed to enforce this.
[opticAxis_t, Rz] = rotateToZeroX(opticAxis_t);

% Definition of alpha consistent with analytic solution
alpha = atan(opticAxis_t(2)/opticAxis_t(3));

nInc = n_inc; % either o or e n parsed above.

% Rotate kinc and the surface normal if alpha was rotated
kinc = Rz*k_inc;
ada = Rz*interaction.normal;

% We know incident medium is isotropic so can define epsilon using
% identity matrix
epsilon = nInc^2*eye(3,3);

% Here assuming uniaxial with alpha paramterizing the crystal axis.
epsilonP = [nO_t.^2, 0, 0 ;  ... 
 0, (nO_i*cos(alpha))^2+(nE_i*sin(alpha))^2, 0.5*(nE_i.^2-nO_i.^2)*sin(2*alpha) ; ...
 0, 0.5*(nE_t.^2-nO_i.^2)*sin(2*alpha), (nE_t*cos(alpha)).^2+(nO_t*sin(alpha)).^2];


% ThetaX/Y in the rotated coordinate system
thetaX = atan(kinc(1)/ kinc(3));
thetaY = atan(kinc(2)/ kinc(3));

% I solved for the analytically using mathematica using 
% PLAOS eqns 19.17 and 19.18.
np = calcNp(thetaY,thetaX,alpha,nO_t,nE_t);


% Transmission waves
k_tO = snell_vector(kinc, ada, nInc, nO_t);
k_tE = snell_vector(kinc, ada, nInc, np);




% Need to compute S

M_tO = calcM(nO_t,nE_t,nO_t,thetaX,thetaY, alpha);

[U,S,V] =svd(M_tO);

% Compute E, H, S using SVD.
E_O = V(:,end);

K_tO = makeK(k_tO);

H_O = nO_t*K_tO*E_O;

S_O = cross(E_O,H_O);
S_O = S_O/norm(S_O);

% Go back to global coordinates now that we have ne
k_tO = inv(Rz)*k_tO;
S_O  = inv(Rz)*S_O;
E_O  = inv(Rz)*E_O;
H_O  = inv(Rz)*H_O;


%% Extraordinary result
% same as ordinary but use np(aka ne) as the index
ne = np;
M_tE = calcM(nO_t,nE_t,np,thetaX,thetaY,alpha);
[U,S,V] =svd(M_tE);

E_E = V(:,end);

K_tE = makeK(k_tE);

H_E = ne*K_tE*E_E;
S_E = cross(E_E,H_E);
S_E = S_E/norm(S_E);

% Go back to global coordinates
k_tE = inv(Rz)*k_tE;
S_E  = inv(Rz)*S_E;
E_E  = inv(Rz)*E_E;
H_E  = inv(Rz)*H_E;



%% Reflection 
% need to account for the possibility of two reflected waves for
% extraordinary wave.  

axisData = mediumIn.AxisData{1};
opticAxisIn = axisData.opticAxis;


a = opticAxisIn(:) / norm(opticAxisIn); % norm here is unnecessary since it is handled when sout.opticAxis is defined
% I need to construct a basis of 3 vectors, 2 orthogonal to a.  It doesn't
% matter which two vectors I choose; the a*a' accomplishes this.
epsilon = nO_i^2*eye(3) + (nE_i^2 - nO_i^2)*(a*a.');

% wave vector parallel to the interface.  TODO:  Update this for arbitrary ada
beta = n_inc * k_inc(1:2);

% Compute o and e reflected waves
[n_ref, k_ref, q_ref] = extraordinaryReflectedQ(beta, opticAxisIn, nO_i, nE_i);
qz_o = -sqrt(nO_i^2 - beta(1)^2 - beta(2)^2);
q_o = [beta(1); beta(2); qz_o];
k_o = q_o / nO_i;

K_re = makeK(k_ref);
M_re = epsilon + (n_ref*K_re)*(n_ref*K_re);
det(M_re) ;% should be close to eps
%norm((epsilon + makeK(q_ref)^2) * E_re) % should be close to eps

[Utmp,Stmp,V] = svd(M_re);
% if opticAxis, ada, and k_inc all in same direction, have a degeneracy
% that needs to be handled.
if all(opticAxis_t == k_inc) && all(opticAxisIn == ada)
   E_re = V(:,end-1);
else
   E_re = V(:,end);
end

H_re = n_ref*K_re*E_re;

S_re = cross(E_re,H_re);
S_re = S_re/norm(S_re);

K_ro = makeK(k_o);
M_ro = epsilon + (nO_i*K_ro)*(nO_i*K_ro);

det(M_ro) ; % should be close to eps

[Utmp,Stmp,V] = svd(M_ro);
%
E_ro = V(:,end);
H_ro = nO_i*K_ro*E_ro;

S_ro = cross(E_ro,H_ro);
S_ro = S_ro/norm(S_ro);



%% -------------------------------------------------------
% s-basis vectors - vectors transerse to the interface
%% -------------------------------------------------------
if (k_inc(1) == 0 && k_inc(2) == 0 && k_inc(3) == 1) % special case normal incidence
    s1 = [1;0;0];
else
    s1 = cross(k_inc, ada);
    s1 = s1/norm(s1);
end
s2 = cross(ada, s1);
s2 = s2/norm(s2);


%% -------------------------------------------------------
% Construct F matrix for this interface
%% -------------------------------------------------------
F = [ ...
 dot(s1,E_ta) dot(s1,E_tb) -dot(s1,E_ro) -dot(s1,E_re);
 dot(s2,E_ta) dot(s2,E_tb) -dot(s2,E_ro) -dot(s2,E_re);
 dot(s1,H_ta) dot(s1,H_tb) -dot(s1,H_ro) -dot(s1,H_re);
 dot(s2,H_ta) dot(s2,H_tb) -dot(s2,H_ro) -dot(s2,H_re) ];


Cm = [ dot(s1,E_inc); dot(s2,E_inc); dot(s1,H_inc); dot(s2,H_inc) ];
% This is the case IV condition from chapter 19
Cn = zeros(size(Cm));


Am = F\Cm; % or F\Cs
Ap = zeros(size(Am)); % or F\Cp 

P_to = [Am(1)*E_ta, zeros(size(E_ta)), S_out]*inv([E_inc, E_n, S_inc]);
P_te = [Am(2)*E_tb, zeros(size(E_ta)), S_out]*inv([E_inc, E_n, S_inc]);

P_ro  = [Am(3)*E_ro, zeros(size(E_ta)), S_ro]*inv([E_inc, E_n, S_inc]);
P_re  = [Am(4)*E_re, zeros(size(E_ta)), S_re]*inv([E_inc, E_n, S_inc]);


% Diagnostic
E_out =   Am(1)*E_ta;
E_out_b = Am(2)*E_tb;
H_out =   Am(1)*H_ta;
H_out_b = Am(2)*H_tb;

Iout_a  = dot(cross(E_out, conj(H_out)), ada);
Iout_b  = dot(cross(E_out_b, conj(H_out_b)), ada);
Iout = Iout_a + Iout_b;
fieldE_out = P_to * ray.fieldE;
fieldScale = prtModalFieldScale(E_out, fieldE_out);
fieldH_out = fieldScale * H_out;
IoutField = dot(cross(fieldE_out, conj(fieldH_out)), ada);

Iout/Iin; % This should not be > 1

% sout_int.P_t = P_t;
% sout_int.E_tb = E_tb; % p 
% sout_int.E_ta = E_ta; % s
% sout_int.k_out = k_out; 


%% Coordinate transformation matrices
[p1,s1] = calcPandSUnitVectors(kinc, interaction.normal);
Oin = calcO(s1,p1,kinc);

[p_o, s_o] = calcPandSUnitVectors(k_tO, interaction.normal);
Oout_o = calcO(s_o,p_o,k_tO);

% For the extraordinary ray, the output local frame follows the energy ray.
[pe, se] = calcPandSUnitVectors(S_E, interaction.normal);
Oout_e = calcO(se,pe,S_E);

Q_o = Oout_o / Oin;
Q_e = Oout_e / Oin;

%% Child rays
childO = makeChildTemplate(ray, hit, normal, mediumOut, "ordinary");
childO.k = k_tO;
childO.S = S_O;
childO.modeE = E_O;
childO.modeH = H_O;
childO.fieldE = P_to * ray.fieldE;
scaleO = prtModalFieldScale(E_O, childO.fieldE);
childO.fieldH = scaleO * H_O;
childO.E = childO.fieldE;
childO.H = childO.fieldH;
childO.P = P_to * ray.P;
childO.Q = Q_o * ray.Q;
childO.O = Oout_o;
childO.localBasis = struct('s', s_o, 'p', p_o, 'basisDirection', k_tO);
childO.amplitude = norm(childO.fieldE);
childO.flux = real(dot(cross(childO.fieldE, conj(childO.fieldH)), interaction.normal));
childO.metadata = struct( ...
    'n', nO_t, ...
    'phaseIndex', nO_t, ...
    'modalScale', scaleO, ...
    'P_interface', P_to, ...
    'P_beforePropagation', ray.P, ...
    'propagationProjector', SD_o, ...
    'epsilon', epsilonP, ...
    'opticAxis', inv(Rz)*opticAxis_t);


childE = makeChildTemplate(ray, hit, normal, mediumOut, "extraordinary");
childE.k = k_tE;
childE.S = S_E;
childE.modeE = E_E;
childE.modeH = H_E;
childE.fieldE = P_te * ray.fieldE;
scaleE = prtModalFieldScale(E_E, childE.fieldE);
childE.fieldH = scaleE * H_E;
childE.E = childE.fieldE;
childE.H = childE.fieldH;
childE.P = P_te * ray.P;
childE.Q = Q_e * ray.Q;
childE.O = Oout_e;
childE.localBasis = struct('s', se, 'p', pe, 'basisDirection', S_E);
childE.amplitude = norm(childE.fieldE);
childE.flux = real(dot(cross(childE.fieldE, conj(childE.fieldH)), interaction.normal));
childE.metadata = struct( ...
    'n', ne, ...
    'n_SE', ne*k_tE'*S_E, ...
    'phaseIndex', ne, ...
    'modalScale', scaleE, ...
    'P_interface', P_te, ...
    'P_beforePropagation', ray.P, ...
    'propagationProjector', SD_e, ...
    'epsilon', epsilonP, ...
    'opticAxis', inv(Rz)*opticAxis_t);

interaction.children = [childO; childE];

%% TODO:  Why is it recalculating E and H?
if surfaceReflective
    reflected_o = makeChildTemplate(ray, hit, normal, mediumIn, "reflected");
    reflected_o.k = k_ro;
    reflected_o.S = S_ro;
    reflected_o.fieldE = P_ro * ray.fieldE;
    reflected_o.fieldH = nr * K_ro * reflected_o.fieldE;
    reflected_o.modeE = reflected_o.fieldE / norm(reflected_o.fieldE);
    reflected_o.modeH = nr * K_ro * reflected_o.modeE;
    reflected_o.E = reflected_o.fieldE;
    reflected_o.H = reflected_o.fieldH;
    reflected_o.P = P_ro * ray.P;
    [pRef, sRef] = interfaceSPBasis(k_r, ada, s1);
    Oref = calcO(sRef, pRef, k_ro);
    reflected_o.Q = (Oref / Oin) * ray.Q;
    reflected_o.O = Oref;
    reflected_o.localBasis = struct('s', sRef, 'p', pRef, 'basisDirection', k_r);
    reflected_o.amplitude = norm(reflected_o.fieldE);
    reflected_o.flux = abs(real(dot(cross(reflected_o.fieldE, conj(reflected_o.fieldH)), ada)));
    reflected_o.metadata = struct( ...
        'n', nr, ...
        'P_interface', P_ro, ...
        'isPropagating', reflected_o.active);

    reflected_e = makeChildTemplate(ray, hit, normal, mediumIn, "reflected");
    reflected_e.k = k_ref;
    reflected_e.S = S_re;
    reflected_e.fieldE = P_re * ray.fieldE;
    reflected_e.fieldH = n_ref * K_re * reflected_e.fieldE;
    reflected_e.modeE = reflected_e.fieldE / norm(reflected_e.fieldE);
    reflected_e.modeH = n_ref * K_re * reflected_e.modeE;
    reflected_e.E = reflected_e.fieldE;
    reflected_e.H = reflected_e.fieldH;
    reflected_e.P = P_re * ray.P;
    [pRef, sRef] = interfaceSPBasis(k_ref, ada, s1);
    Oref = calcO(sRef, pRef, k_ref);
    reflected_e.Q = (Oref / Oin) * ray.Q;
    reflected_e.O = Oref;
    reflected_e.localBasis = struct('s', sRef, 'p', pRef, 'basisDirection', k_r);
    reflected_e.amplitude = norm(reflected_e.fieldE);
    reflected_e.flux = abs(real(dot(cross(reflected_e.fieldE, conj(reflected_e.fieldH)), ada)));
    reflected_e.metadata = struct( ...
        'n', n_ref, ...
        'P_interface', P_re, ...
        'isPropagating', reflected_e.active);    
    children = [children; reflected_e];
end

end

function M = calcM(no,nE, n,tX,tY,alpha)

% helper expressions
Tx = tan(tX);
Ty = tan(tY);
A  = abs(1 + Tx^2 + Ty^2);        % original Abs(...)
Tsum = Tx^2 + Ty^2;
commonSqrt = sqrt(n^2*A - Tsum); % sqrt(n^2*A - (tan^2(tX)+tan^2(tY)))

% simplified matrix
M = [ (no^2 - n^2) + Tx^2 / A, ...
      (Tx*Ty) / A, ...
      Tx * commonSqrt / A;
      
      (Tx*Ty) / A, ...
      (no^2*cos(alpha)^2 + nE^2*sin(alpha)^2 - n^2) + Ty^2 / A, ...
      0.5*(nE^2 - no^2)*sin(2*alpha) + Ty * commonSqrt / A;
      
      Tx * commonSqrt / A, ...
      0.5*(nE^2 - no^2)*sin(2*alpha) + Ty * commonSqrt / A, ...
      nE^2*cos(alpha)^2 + no^2*sin(alpha)^2 - Tsum / A ];

end

function nprime = calcNp(tY,tX,alpha,no,nE)

n1 = ( sqrt( ...
    21*nE^4 - 30*nE^2*no^2 + 9*no^4 + ...
    12*nE^4*cos(2*tY) - 8*nE^2*no^2*cos(2*tY) - 4*no^4*cos(2*tY) - ...
    9*nE^4*cos(4*tY) + 22*nE^2*no^2*cos(4*tY) - 13*no^4*cos(4*tY) + ...
    24*nE^4*cos(2*alpha) - 12*nE^2*no^2*cos(2*alpha) - 12*no^4*cos(2*alpha) + ...
    16*nE^4*cos(2*tY)*cos(2*alpha) - 16*nE^2*no^2*cos(2*tY)*cos(2*alpha) - ...
    8*nE^4*cos(4*tY)*cos(2*alpha) - 4*nE^2*no^2*cos(4*tY)*cos(2*alpha) + ...
    12*no^4*cos(4*tY)*cos(2*alpha) + ...
    256*nE^2*no^2*abs(sec(tY)^2 + tan(tX)^2).*cos(tX).^4.*cos(tY).^4.* ...
        (nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)) + ...
    3*nE^4*cos(4*alpha) - 6*nE^2*no^2*cos(4*alpha) + 3*no^4*cos(4*alpha) + ...
    4*nE^4*cos(2*tY)*cos(4*alpha) - 8*nE^2*no^2*cos(2*tY)*cos(4*alpha) + ...
    4*no^4*cos(2*tY)*cos(4*alpha) + ...
    nE^4*cos(4*tY)*cos(4*alpha) - 2*nE^2*no^2*cos(4*tY)*cos(4*alpha) + ...
    no^4*cos(4*tY)*cos(4*alpha) - ...
    2*(nE^2 - no^2)*cos(4*tX).*cos(tY).^2 .* ( ...
        -2*nE^2 + 10*no^2 + ...
        2*(7*nE^2 - 3*no^2)*cos(2*tY) + (nE^2 - no^2)*cos(2*tY - 4*alpha) + ...
        8*nE^2*cos(2*tY - 2*alpha) + 4*no^2*cos(2*tY - 2*alpha) - ...
        8*no^2*cos(2*alpha) + 2*nE^2*cos(4*alpha) - 2*no^2*cos(4*alpha) + ...
        8*nE^2*cos(2*(tY + alpha)) + 4*no^2*cos(2*(tY + alpha)) + ...
        nE^2*cos(2*tY + 4*alpha) - no^2*cos(2*tY + 4*alpha)) + ...
    32*nE^4*cos(2*tX).*sin(2*tY).^2 - 64*nE^2*no^2*cos(2*tX).*sin(2*tY).^2 + ...
    32*no^4*cos(2*tX).*sin(2*tY).^2 + ...
    32*nE^4*cos(2*tX).*cos(2*alpha).*sin(2*tY).^2 - ...
    32*no^4*cos(2*tX).*cos(2*alpha).*sin(2*tY).^2 + ...
    128*sqrt(2)*sqrt( ...
        no^2*(nE^2 - no^2)^2 .* cos(tX).^6 .* cos(tY).^4 .* sin(tY).^2 .* ( ...
            -3*nE^2 - no^2 - nE^2*cos(2*alpha) + no^2*cos(2*alpha) + ...
            4*nE^2*abs(sec(tY)^2 + tan(tX)^2).*cos(tX).^2.*cos(tY).^2 .* ...
                (nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)) + ...
            2*nE^2*cos(2*tY).*sin(alpha).^2 - 2*no^2*cos(2*tY).*sin(alpha).^2 + ...
            cos(2*tX) .* ( cos(2*tY)*(3*nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)) + ...
                           2*(-nE^2 + no^2)*sin(alpha).^2 ) ...
        ) .* sin(2*alpha).^2 ) ...
)) ./ ( 8*sqrt(2)*sqrt(abs(sec(tY)^2 + tan(tX)^2)) .* ...
        sqrt( cos(tX).^4 .* cos(tY).^4 .* (nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)).^2 ) );

n2 = ( sqrt( ...
    21*nE^4 - 30*nE^2*no^2 + 9*no^4 + ...
    12*nE^4*cos(2*tY) - 8*nE^2*no^2*cos(2*tY) - 4*no^4*cos(2*tY) - ...
    9*nE^4*cos(4*tY) + 22*nE^2*no^2*cos(4*tY) - 13*no^4*cos(4*tY) + ...
    24*nE^4*cos(2*alpha) - 12*nE^2*no^2*cos(2*alpha) - 12*no^4*cos(2*alpha) + ...
    16*nE^4*cos(2*tY)*cos(2*alpha) - 16*nE^2*no^2*cos(2*tY)*cos(2*alpha) - ...
    8*nE^4*cos(4*tY)*cos(2*alpha) - 4*nE^2*no^2*cos(4*tY)*cos(2*alpha) + ...
    12*no^4*cos(4*tY)*cos(2*alpha) + ...
    256*nE^2*no^2*abs(sec(tY)^2 + tan(tX)^2).*cos(tX).^4.*cos(tY).^4.* ...
        (nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)) + ...
    3*nE^4*cos(4*alpha) - 6*nE^2*no^2*cos(4*alpha) + 3*no^4*cos(4*alpha) + ...
    4*nE^4*cos(2*tY)*cos(4*alpha) - 8*nE^2*no^2*cos(2*tY)*cos(4*alpha) + ...
    4*no^4*cos(2*tY)*cos(4*alpha) + ...
    nE^4*cos(4*tY)*cos(4*alpha) - 2*nE^2*no^2*cos(4*tY)*cos(4*alpha) + ...
    no^4*cos(4*tY)*cos(4*alpha) - ...
    2*(nE^2 - no^2).*cos(4*tX).*cos(tY).^2 .* ( ...
        -2*nE^2 + 10*no^2 + ...
        2*(7*nE^2 - 3*no^2).*cos(2*tY) + (nE^2 - no^2).*cos(2*tY - 4*alpha) + ...
        8*nE^2*cos(2*tY - 2*alpha) + 4*no^2*cos(2*tY - 2*alpha) - ...
        8*no^2*cos(2*alpha) + 2*nE^2*cos(4*alpha) - 2*no^2*cos(4*alpha) + ...
        8*nE^2*cos(2*(tY + alpha)) + 4*no^2*cos(2*(tY + alpha)) + ...
        nE^2*cos(2*tY + 4*alpha) - no^2*cos(2*tY + 4*alpha) ) + ...
    32*nE^4*cos(2*tX).*sin(2*tY).^2 - 64*nE^2*no^2*cos(2*tX).*sin(2*tY).^2 + ...
    32*no^4*cos(2*tX).*sin(2*tY).^2 + ...
    32*nE^4*cos(2*tX).*cos(2*alpha).*sin(2*tY).^2 - ...
    32*no^4*cos(2*tX).*cos(2*alpha).*sin(2*tY).^2 - ...
    128*sqrt(2)*sqrt( ...
        no^2*(nE^2 - no^2).^2 .* cos(tX).^6 .* cos(tY).^4 .* sin(tY).^2 .* ( ...
            -3*nE^2 - no^2 - nE^2*cos(2*alpha) + no^2*cos(2*alpha) + ...
            4*nE^2*abs(sec(tY)^2 + tan(tX)^2).*cos(tX).^2.*cos(tY).^2 .* ...
                (nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)) + ...
            2*nE^2*cos(2*tY).*sin(alpha).^2 - 2*no^2*cos(2*tY).*sin(alpha).^2 + ...
            cos(2*tX).*( cos(2*tY).*(3*nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)) + ...
                         2*(-nE^2 + no^2).*sin(alpha).^2 ) ...
        ) .* sin(2*alpha).^2 ) ...
)) ./ ( 8*sqrt(2).*sqrt(abs(sec(tY)^2 + tan(tX)^2)) .* ...
        sqrt( cos(tX).^4 .* cos(tY).^4 .* (nE^2 + no^2 + (nE^2 - no^2)*cos(2*alpha)).^2 ) );

% Compute angle between tX,tY and alpha
k_alpha = [0; sin(alpha); cos(alpha)];
% Following convention in chapter 11 of PLAOS
kincNorm = abs(sqrt(tan(tX)^2+tan(tY)^2+1));
kinc = 1/kincNorm*[tan(tX); tan(tY); 1];

theta = acos(kinc(3));

% This is to check if the two vectors are in the same quadrant or not in the yz plane.

% First, figure out which solution is closest to the extraordinatry index
n1tst = abs(n1-nE);
n2tst = abs(n2-nE);

if n1tst < n2tst
    ncE = n1;
    ncO = n2;
else
    ncE = n2;
    ncO = n1;
end
% if (dot(k_alpha,kinc)) < 0
%     nprime = ncE; % if positive crystal axis then nprime should be closer to nO
% else
%     nprime = ncO;
% end

if (sign(k_alpha(2)*kinc(2)) < 0) || (sign(k_alpha(3)*kinc(3)) < 0)
    %nprime = n2;
    nprime = ncE; % if positive crystal axis then nprime should be closer to nO
else
    %nprime = n1;
    nprime = ncO; % if negative crystal axis then nprime should be closer to nE
end

% if (sign(k_alpha(2)*kinc(2)) < 0) || (sign(k_alpha(3)*kinc(3)) > 0)
%     %nprime = n2;
% else
%     %nprime = n1;
% end



end

