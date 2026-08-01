% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function interaction = traceUniaxialToUniaxial(interaction, mediumIn, mediumOut, ray, hit, normal, options)
%TRACEUNIAXIALTOUNIAXIAL trace rays from uniaxial to uniaxial materials
%Based on PLAOS Section 19.7.4


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



% Conserve the tangential wavevector and solve the two output-medium
% dielectric eigenconditions directly. This remains valid for arbitrary
% interface normals and optic-axis orientations.
qInc = n_inc*k_inc;
qTangential = qInc - dot(qInc, ada)*ada;
[modeO, modeE] = uniaxialModesFromTangentialQ( ...
    qTangential, ada, nO_t, nE_t, opticAxis_t, 1);
k_tO = modeO.k;
S_O = modeO.S;
E_O = modeO.E;
H_O = modeO.H;
k_tE = modeE.k;
S_E = modeE.S;
E_E = modeE.E;
H_E = modeE.H;
ne = modeE.phaseIndex;
epsilonP = nO_t^2*eye(3) + ...
    (nE_t^2 - nO_t^2)*(opticAxis_t*opticAxis_t.');


%% Reflection 
% need to account for the possibility of two reflected waves for
% extraordinary wave.  

axisData = mediumIn.AxisData{1};
opticAxisIn = axisData.opticAxis;


a = opticAxisIn(:) / norm(opticAxisIn); % norm here is unnecessary since it is handled when sout.opticAxis is defined
% I need to construct a basis of 3 vectors, 2 orthogonal to a.  It doesn't
% matter which two vectors I choose; the a*a' accomplishes this.
epsilon = nO_i^2*eye(3) + (nE_i^2 - nO_i^2)*(a*a.');

% Compute o and e reflected waves while conserving tangential q.
[n_ref, k_ref, q_ref] = extraordinaryReflectedQ( ...
    qInc, ada, opticAxisIn, nO_i, nE_i);
qParallel = qInc - dot(qInc,ada)*ada;
qNormalO = -sqrt(nO_i^2 - dot(qParallel,qParallel));
qRefO = qParallel + qNormalO*ada;
k_ro = qRefO / nO_i;

K_re = makeK(k_ref);
M_re = epsilon + (n_ref*K_re)*(n_ref*K_re);
det(M_re) ;% should be close to eps
%norm((epsilon + makeK(q_ref)^2) * E_re) % should be close to eps

[Utmp,Stmp,V] = svd(M_re);
% if opticAxis, ada, and k_inc all in same direction, have a degeneracy
% that needs to be handled.
if all(opticAxisIn == k_inc) && all(opticAxisIn == ada)
   E_re = V(:,end-1);
else
   E_re = V(:,end);
end

H_re = n_ref*K_re*E_re;

S_re = cross(E_re,H_re);
S_re = S_re/norm(S_re);

K_ro = makeK(k_ro);
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
 dot(s1,E_O) dot(s1,E_E) -dot(s1,E_ro) -dot(s1,E_re);
 dot(s2,E_O) dot(s2,E_E) -dot(s2,E_ro) -dot(s2,E_re);
 dot(s1,H_O) dot(s1,H_E) -dot(s1,H_ro) -dot(s1,H_re);
 dot(s2,H_O) dot(s2,H_E) -dot(s2,H_ro) -dot(s2,H_re) ];

Cm = [ dot(s1,E_inc); dot(s2,E_inc); dot(s1,H_inc); dot(s2,H_inc) ];
% This is the case IV condition from chapter 19
Cn = zeros(size(Cm));


Am = F\Cm; % or F\Cs
Ap = zeros(size(Am)); % or F\Cp 

P_to = [Am(1)*E_O, zeros(size(E_O)), S_O]*inv([E_inc, E_n, S_inc]);
P_te = [Am(2)*E_E, zeros(size(E_O)), S_E]*inv([E_inc, E_n, S_inc]);

P_ro  = [Am(3)*E_ro, zeros(size(E_O)), S_ro]*inv([E_inc, E_n, S_inc]);
P_re  = [Am(4)*E_re, zeros(size(E_O)), S_re]*inv([E_inc, E_n, S_inc]);

SD_o = k_tO * transpose(k_inc);
SD_e = S_E * transpose(k_inc);

%% Diagnostics

% First let's check flux ratios.  If there is no absorption we should be
% able to prove conversation of energy if everything is set up correctly.
% If we have E_s input, then have to look at the intensity for transmitted
% o and e and reflected s and p.  

Eout_o = Am(1)*E_O;
Hout_o = Am(1)*H_O;

Eout_e = Am(2)*E_E;
Hout_e = Am(2)*H_E;

Eref_o = Am(3)*E_ro;
Href_o = Am(3)*H_ro;

Eref_e = Am(4)*E_re;
Href_e = Am(4)*H_re;

Iout_o = dot(cross(Eout_o, conj(Hout_o)), ada);
Iout_e = dot(cross(Eout_e, conj(Hout_e)), ada);

Iout_ro = dot(cross(Eref_o, conj(Href_o)), -ada);
Iout_re = dot(cross(Eref_e, conj(Href_e)), -ada);


Iout_all = Iout_o+Iout_e + Iout_ro+Iout_re;


transmissionFluxRatio = Iout_all/Iin; % This should not be == 1 if there is no absorption

% We can also detect direction vectors. 
% E dot H should be 0
% Check transverse projections of E are equal across interface

% Diagnostic
E_out =   Am(1)*E_O;
E_out_b = Am(2)*E_E;
H_out =   Am(1)*H_O;
H_out_b = Am(2)*H_E;

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
[p1,s1] = calcPandSUnitVectors(k_inc, interaction.normal);
Oin = calcO(s1,p1,k_inc);

[p_o, s_o] = calcPandSUnitVectors(k_tO, interaction.normal);
Oout_o = calcO(s_o,p_o,k_tO);

% For the extraordinary ray, the output local frame follows the energy ray.
[pe, se] = calcPandSUnitVectors(S_E, interaction.normal);
Oout_e = calcO(se,pe,S_E);

Q_o = Oout_o / Oin;
Q_e = Oout_e / Oin;

%% Child rays
childO = buildPolarizedChild(ray, hit, normal, mediumOut, ...
    Mode="ordinary", BranchType="transmitted", ...
    k=k_tO, S=S_O, P=P_to, Q=Q_o, O=Oout_o, ...
    LocalBasis=struct('s', s_o, 'p', p_o, 'basisDirection', k_tO), ...
    Index=nO_t, ModeE=E_O, ModeH=H_O, ...
    FluxNormal=interaction.normal, PropagationProjector=SD_o, ...
    Metadata=struct( ...
    'epsilon', epsilonP, ...
    'opticAxis', opticAxis_t));


childE = buildPolarizedChild(ray, hit, normal, mediumOut, ...
    Mode="extraordinary", BranchType="transmitted", ...
    k=k_tE, S=S_E, P=P_te, Q=Q_e, O=Oout_e, ...
    LocalBasis=struct('s', se, 'p', pe, 'basisDirection', S_E), ...
    Index=ne, ModeE=E_E, ModeH=H_E, ...
    FluxNormal=interaction.normal, PhaseIndex=ne, ...
    PropagationProjector=SD_e, Metadata=struct( ...
    'n_SE', ne*k_tE'*S_E, ...
    'epsilon', epsilonP, ...
    'opticAxis', opticAxis_t));

interaction.children = [childO; childE];

[pRefO, sRefO] = calcPandSUnitVectors(k_ro, ada, s1);
OrefO = calcO(sRefO, pRefO, k_ro);
QrefO = OrefO / Oin;
SD_ro = S_ro * transpose(k_inc);
reflected_o = buildPolarizedChild(ray, hit, normal, mediumIn, ...
    Mode="ordinary", BranchType="reflected", ...
    k=k_ro, S=S_ro, P=P_ro, Q=QrefO, O=OrefO, ...
    LocalBasis=struct('s', sRefO, 'p', pRefO, 'basisDirection', k_ro), ...
    Index=nO_i, ModeE=E_ro, ModeH=H_ro, FluxNormal=-ada, ...
    PropagationProjector=SD_ro, Metadata=struct( ...
    'epsilon', epsilon, 'opticAxis', opticAxisIn));

[pRefE, sRefE] = calcPandSUnitVectors(k_ref, ada, s1);
OrefE = calcO(sRefE, pRefE, k_ref);
QrefE = OrefE / Oin;
SD_re = S_re * transpose(k_inc);
reflected_e = buildPolarizedChild(ray, hit, normal, mediumIn, ...
    Mode="extraordinary", BranchType="reflected", ...
    k=k_ref, S=S_re, P=P_re, Q=QrefE, O=OrefE, ...
    LocalBasis=struct('s', sRefE, 'p', pRefE, 'basisDirection', k_ref), ...
    Index=n_ref, ModeE=E_re, ModeH=H_re, FluxNormal=-ada, ...
    PhaseIndex=n_ref, PropagationProjector=SD_re, Metadata=struct( ...
    'n_SE', n_ref*dot(k_ref,S_re), ...
    'epsilon', epsilon, 'opticAxis', opticAxisIn));
interaction.children = ...
    [interaction.children; reflected_o; reflected_e];

interaction.frames.Oin = Oin;
interaction.frames.OoutOrdinary = Oout_o;
interaction.frames.OoutExtraordinary = Oout_e;
interaction.frames.OrefOrdinary = OrefO;
interaction.frames.OrefExtraordinary = OrefE;
interaction.P = struct( ...
    'transmittedOrdinary', P_to, ...
    'transmittedExtraordinary', P_te, ...
    'reflectedOrdinary', P_ro, ...
    'reflectedExtraordinary', P_re);
interaction.Q = struct( ...
    'transmittedOrdinary', Q_o, ...
    'transmittedExtraordinary', Q_e, ...
    'reflectedOrdinary', QrefO, ...
    'reflectedExtraordinary', QrefE);
interaction.coefficients = struct( ...
    'Am', Am, ...
    'Ap', Ap, ...
    'transmittedOrdinary', Am(1), ...
    'transmittedExtraordinary', Am(2), ...
    'reflectedOrdinary', Am(3), ...
    'reflectedExtraordinary', Am(4));
interaction.diagnostics = struct( ...
    'F', F, ...
    'Cm', Cm, ...
    'Cn', Cn, ...
    'boundaryResidual_m', F*Am - Cm, ...
    'boundaryResidual_n', F*Ap - Cn, ...
    'Iin', Iin, ...
    'IinField', IinField, ...
    'totalFluxRatio', transmissionFluxRatio, ...
    'totalFluxRatioField', ...
        sum([interaction.children.flux])/IinField);

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
