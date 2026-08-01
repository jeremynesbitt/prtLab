% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function interaction = traceIsotropicToUniaxial(interaction, mediumIn, mediumOut, ray, hit, normal, options)
%TRACEISOTROPICTOUNIAXIAL Placeholder for isotropic-to-uniaxial tracing.
% Jeremy Nesbitt
% October 2025
% Rev 2 - fix bug with calculating wrong nprime for negative alpha
% and use direction cosines as input
% Rev 3 - add reflection coefficients
% Rev 4 - incorporate coordinate transformation matrices
% Rev 5 - update for prtLab interface

axisData = mediumOut.AxisData{1};
opticAxis = axisData.opticAxis;
if norm(opticAxis) ~= 1
    opticAxis = opticAxis/norm(opticAxis);
end


deg = pi/180;

kinc = interaction.incident.k;

% The analytic solutions below require alpha to be in the yz plane.
% rotate the coordinate system if needed to enforce this.
[opticAxis, Rz] = rotateToZeroX(opticAxis);

% Definition of alpha consistent with analytic solution
alpha = atan(opticAxis(2)/opticAxis(3));

iStruct = mediumIn.IndexData{1};
oStruct = mediumOut.IndexData{1};

nInc = iStruct.n;

nO = oStruct.nO;
nE = oStruct.nE;


% Rotate kinc and the surface normal if alpha was rotated
kinc = Rz*kinc;
ada = Rz*interaction.normal;

% We know incident medium is isotropic so can define epsilon using
% identity matrix
epsilon = nInc^2*eye(3,3);

% Here assuming uniaxial with alpha paramterizing the crystal axis.
epsilonP = [nO.^2, 0, 0 ;  ... 
 0, (nO*cos(alpha))^2+(nE*sin(alpha))^2, 0.5*(nE.^2-nO.^2)*sin(2*alpha) ; ...
 0, 0.5*(nE.^2-nO.^2)*sin(2*alpha), (nE*cos(alpha)).^2+(nO*sin(alpha)).^2];


% ThetaX/Y in the rotated coordinate system
thetaX = atan(kinc(1)/ kinc(3));
thetaY = atan(kinc(2)/ kinc(3));

% I solved for the analytically using mathematica using 
% PLAOS eqns 19.17 and 19.18.
np = calcNp(thetaY,thetaX,alpha,nO,nE);


% Transmission waves
k_tO = snell_vector(kinc, ada, nInc, nO);
k_tE = snell_vector(kinc, ada, nInc, np);




% Need to compute S

M_tO = calcM(nO,nE,nO,thetaX,thetaY, alpha);

[U,S,V] =svd(M_tO);

% Compute E, H, S using SVD.
E_O = V(:,end);

K_tO = makeK(k_tO);

H_O = nO*K_tO*E_O;

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
M_tE = calcM(nO,nE,np,thetaX,thetaY,alpha);
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

nr = nInc;
k_r = kinc-2*(kinc'*ada)*ada;
K_r = makeK(k_r);

M_r = epsilon + K_r*K_r;

[U,S,V] =svd(M_r);

Er_s = V(:,end);
Er_p = -V(:,end-1); % Again having sign error issues TODO:  Close on this
Hr_s = nr*K_r*Er_s;
Hr_p = nr*K_r*Er_p;

% Doesn't matter whether I choose s or p here.  they point in the same
% direction
Sr = cross(Er_p,Hr_p);
Sr = Sr/norm(Sr);


% Go back to global coordinates
k_r = inv(Rz)*k_r;
Sr  = inv(Rz)*Sr;
Er_s  = inv(Rz)*Er_s;
Er_p  = inv(Rz)*Er_p;
Hr_s  = inv(Rz)*Hr_s;
Hr_p  = inv(Rz)*Hr_p;
K_r   = inv(Rz)*K_r;

% revert kinc since we are done with the rotated coordinate system.
kinc = inv(Rz)*kinc;

%% P matrices
E_to = E_O;
E_te = E_E;
H_to = H_O;
H_te = H_E;
E_rs = Er_s;
E_rp = Er_p;
H_rs = Hr_s;
H_rp = Hr_p;
S_inc = kinc;

if (kinc(1) == 0 && kinc(2) == 0 && kinc(3) == 1) % special case normal incidence
    s1 = [1;0;0];
else
    s1 = cross(kinc,interaction.normal); 
    s1 = s1/norm(s1);
end
s2 = cross(interaction.normal,s1);
s2 = s2/norm(s2);

%% -------------------------------------------------------
% Construct F matrix (Eq. 19.47)
%% -------------------------------------------------------
F = [ ...
 dot(s1,E_to) dot(s1,E_te) -dot(s1,E_rs) -dot(s1,E_rp);
 dot(s2,E_to) dot(s2,E_te) -dot(s2,E_rs) -dot(s2,E_rp);
 dot(s1,H_to) dot(s1,H_te) -dot(s1,H_rs) -dot(s1,H_rp);
 dot(s2,H_to) dot(s2,H_te) -dot(s2,H_rs) -dot(s2,H_rp) ];


% -------------------------------------------------------
% Incident s and p states
% -------------------------------------------------------

if (kinc(1) == 0 && kinc(2) == 0 && kinc(3) == 1) % special case normal incidence
    Einc_s = [1;0;0];
else
    Einc_s = cross(S_inc,interaction.normal)/norm(cross(S_inc,interaction.normal));
end

Einc_p = cross(S_inc,Einc_s)/norm(cross(S_inc,Einc_s));
Kinc = makeK(kinc);
Hinc_s = nInc * Kinc * Einc_s;
Hinc_p = nInc * Kinc * Einc_p;


Cs = [ dot(s1,Einc_s); dot(s2,Einc_s); dot(s1,Hinc_s); dot(s2,Hinc_s) ];
Cp = [ dot(s1,Einc_p); dot(s2,Einc_p); dot(s1,Hinc_p); dot(s2,Hinc_p) ];


As = F\Cs; % or F\Cs
Ap = F\Cp; % or F\Cp 

%% -------------------------------------------------------
% Amplitude coefficients
%% -------------------------------------------------------
a_s_to = As(1); a_s_te = As(2);
a_s_rs = As(3); a_s_rp = As(4);

a_p_to = Ap(1); a_p_te = Ap(2);
a_p_rs = Ap(3); a_p_rp = Ap(4);

S_to = real(cross(E_to,conj(H_to)))/norm(cross(E_to,conj(H_to)));
S_te = real(cross(E_te,conj(H_te))/norm(cross(E_te,conj(H_te))));

S_rs = real(cross(E_rs,conj(H_rs))/norm(cross(E_rs,conj(H_rs))));

%% -------------------------------------------------------
% Construct 3×3 P matrices (Section 19.7.2)
%% -------------------------------------------------------
P_to = [ a_s_to*E_to, a_p_to*E_to, S_to ]*inv([Einc_s, Einc_p, S_inc]);
P_te = [ a_s_te*E_te, a_p_te*E_te, S_te ]*inv([Einc_s, Einc_p, S_inc]);

P_r  = [ ...
    a_s_rs*E_rs + a_s_rp*E_rp, ...
    a_p_rs*E_rs + a_p_rp*E_rp, ...
    S_rs] / [Einc_s, Einc_p, S_inc];

SD_o = k_tO * transpose(kinc);
SD_e = S_E * transpose(kinc);



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
[p_r, s_r] = calcPandSUnitVectors(k_r, interaction.normal, s1);
Oref = calcO(s_r,p_r,k_r);
Q_r = Oref*(getIReflect / Oin);
SD_r = Sr * transpose(kinc);

%% Diagnostics

% First let's check flux ratios.  If there is no absorption we should be
% able to prove conversation of energy if everything is set up correctly.
% If we have E_s input, then have to look at the intensity for transmitted
% o and e and reflected s and p.  

Eout_so = As(1)*E_O;
Hout_so = As(1)*H_O;

Eout_po = Ap(1)*E_O;
Hout_po = Ap(1)*H_O;

Eout_se = As(2)*E_E;
Hout_se = As(2)*H_E;

Eout_pe = Ap(2)*E_E;
Hout_pe = Ap(2)*H_E;

Eref_ss = As(3)*Er_s;
Href_ss = As(3)*Hr_s;

Eref_ps = Ap(3)*Er_s;
Href_ps = Ap(3)*Hr_s;

Eref_sp = As(4)*Er_p;
Href_sp = As(4)*Hr_p;

Eref_pp = Ap(4)*Er_p;
Href_pp = Ap(4)*Hr_p;

Iin_s = dot(cross(Einc_s, conj(Hinc_s)), ada);
Iin_p = dot(cross(Einc_p, conj(Hinc_p)), ada);

Iout_so = dot(cross(Eout_so, conj(Hout_so)), ada);
Iout_se = dot(cross(Eout_se, conj(Hout_se)), ada);

Iout_po = dot(cross(Eout_po, conj(Hout_po)), ada);
Iout_pe = dot(cross(Eout_pe, conj(Hout_pe)), ada);

Iout_rss = dot(cross(Eref_ss, conj(Href_ss)), -ada);
Iout_rsp = dot(cross(Eref_sp, conj(Href_sp)), -ada);

Iout_rps = dot(cross(Eref_ps, conj(Href_ps)), -ada);
Iout_rpp = dot(cross(Eref_pp, conj(Href_pp)), -ada);


Iout_all_s = Iout_so+Iout_se + Iout_rss+Iout_rsp;

Iout_all_p = Iout_po+Iout_pe + Iout_rps+Iout_rpp;


transmissionFluxRatio_s = Iout_all_s/Iin_s; % This should not be == 1 if there is no absorption

transmissionFluxRatio_p = Iout_all_p/Iin_p; % This should not be == 1 if there is no absorption

% We can also detect direction vectors. 
% E dot H should be 0
% Check transverse projections of E are equal across interface


%% Child rays
childO = buildPolarizedChild(ray, hit, normal, mediumOut, ...
    Mode="ordinary", BranchType="transmitted", ...
    k=k_tO, S=S_O, P=P_to, Q=Q_o, O=Oout_o, ...
    LocalBasis=struct('s', s_o, 'p', p_o, 'basisDirection', k_tO), ...
    Index=nO, ModeE=E_O, ModeH=H_O, ...
    FluxNormal=interaction.normal, PropagationProjector=SD_o, ...
    Metadata=struct( ...
    'epsilon', epsilonP, ...
    'opticAxis', inv(Rz)*opticAxis));

childE = buildPolarizedChild(ray, hit, normal, mediumOut, ...
    Mode="extraordinary", BranchType="transmitted", ...
    k=k_tE, S=S_E, P=P_te, Q=Q_e, O=Oout_e, ...
    LocalBasis=struct('s', se, 'p', pe, 'basisDirection', S_E), ...
    Index=ne, ModeE=E_E, ModeH=H_E, ...
    FluxNormal=interaction.normal, PhaseIndex=ne, ...
    PropagationProjector=SD_e, Metadata=struct( ...
    'n_SE', ne*k_tE'*S_E, ...
    'epsilon', epsilonP, ...
    'opticAxis', inv(Rz)*opticAxis));

childR = buildPolarizedChild(ray, hit, normal, mediumIn, ...
    Mode="isotropic", BranchType="reflected", ...
    k=k_r, S=Sr, P=P_r, Q=Q_r, O=Oref, ...
    LocalBasis=struct('s', s_r, 'p', p_r, 'basisDirection', k_r), ...
    Index=nInc, FluxNormal=-interaction.normal, ...
    PropagationProjector=SD_r);

interaction.children = [childO; childE; childR];

%% Surface-level records
interaction.frames.Oin = Oin;
interaction.frames.Oout_o = Oout_o;
interaction.frames.Oout_e = Oout_e;
interaction.frames.Q_o = Q_o;
interaction.frames.Q_e = Q_e;
interaction.frames.Oref = Oref;
interaction.frames.Qref = Q_r;
interaction.frames.inputBasis = struct('s', s1, 'p', p1, 'basisDirection', kinc);
interaction.frames.outputBasis_o = childO.localBasis;
interaction.frames.outputBasis_e = childE.localBasis;
interaction.frames.reflectedBasis = childR.localBasis;

interaction.P = struct( ...
    'to_o', P_to, ...
    'to_e', P_te, ...
    'reflected', P_r, ...
    'SD_o', SD_o, ...
    'SD_e', SD_e);
interaction.Q = struct( ...
    'ordinary', Q_o, ...
    'extraordinary', Q_e, ...
    'reflected', Q_r);
interaction.coefficients = struct( ...
    'As', As, ...
    'Ap', Ap, ...
    'a_s_to', a_s_to, ...
    'a_s_te', a_s_te, ...
    'a_s_rs', a_s_rs, ...
    'a_s_rp', a_s_rp, ...
    'a_p_to', a_p_to, ...
    'a_p_te', a_p_te, ...
    'a_p_rs', a_p_rs, ...
    'a_p_rp', a_p_rp);
interaction.diagnostics = struct( ...
    'epsilon', epsilonP, ...
    'opticAxis', inv(Rz)*opticAxis, ...
    'F', F, ...
    'Cs', Cs, ...
    'Cp', Cp, ...
    'boundaryResidual_s', F*As - Cs, ...
    'boundaryResidual_p', F*Ap - Cp, ...
    'k_reflected', k_r, ...
    'S_reflected', Sr, ...
    'Er_s', Er_s, ...
    'Er_p', Er_p, ...
    'Hr_s', Hr_s, ...
    'Hr_p', Hr_p);


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
