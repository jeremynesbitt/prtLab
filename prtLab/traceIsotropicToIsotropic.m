% Copyright (c) 2026 Jeremy Nesbitt. All rights reserved.
% Use of this source code is governed by a BSD-style license that can be
% found in the LICENSE file.

function interaction = traceIsotropicToIsotropic(interaction, mediumIn, mediumOut, ray, hit, normal, ~)
%TRACEISOTROPICTOISOTROPIC Trace an isotropic-to-isotropic interface.

% arbitrary interface!
%
% no thickness here yet so not properly propagating to the next surface

%% -------------------------------------------------------
% Geometry and indices
%% -------------------------------------------------------
k_inc = interaction.incident.k;
S_inc = interaction.incident.S;
n1 = incidentIndexFromRay(ray, mediumIn);
n2 = mediumOut.IndexData{1}.n;
ada = interaction.normal;
tir = isTotalInternalReflection(k_inc, ada, n1, n2);

% Use snell's law to calculate the exit fields.


k_out = snell_vector(k_inc, ada, n1, n2);
S_out = k_out;

[Eout_p, Eout_s] = calcPandSUnitVectors(k_out, ada, s1FromRay(k_inc, ada));
K_out = makeK(k_out);
Hout_s = n2 * K_out * Eout_s;
Hout_p = n2 * K_out * Eout_p;


%% Reflection 

nr = n1; 
k_r = k_inc-2*(k_inc'*ada)*ada;

[Er_p, Er_s] = calcPandSUnitVectors(k_r, ada, s1FromRay(k_inc, ada));

K_r = makeK(k_r);
Hr_s = nr*K_r*Er_s;
Hr_p = nr*K_r*Er_p;

% Doesn't matter whether I choose s or p here.  they point in the same
% direction
S_r = cross(Er_p,Hr_p);
S_r = prtNorm(S_r);


%% -------------------------------------------------------
% s-basis vectors - vectors transerse to the interface
%% -------------------------------------------------------
if norm(cross(k_inc, ada)) < 100*eps(max(1, norm(k_inc)*norm(ada))) % special case normal incidence
    s1 = [1;0;0];
else
    s1 = cross(k_inc, ada);
    s1 = s1/norm(s1);
end
s2 = cross(ada, s1);
s2 = s2/norm(s2);

%% -------------------------------------------------------
% Incident s and p states at this interface
%% -------------------------------------------------------
[Einc_p, Einc_s] = calcPandSUnitVectors(k_inc, ada, s1);
K_inc = makeK(k_inc);
Hinc_s = n1 * K_inc * Einc_s;
Hinc_p = n1 * K_inc * Einc_p;

%% -------------------------------------------------------
% Construct F matrix for this interface
%% -------------------------------------------------------
F = [ ...
 dot(s1,Eout_s) dot(s1,Eout_p) -dot(s1,Er_s) -dot(s1,Er_p);
 dot(s2,Eout_s) dot(s2,Eout_p) -dot(s2,Er_s) -dot(s2,Er_p);
 dot(s1,Hout_s) dot(s1,Hout_p) -dot(s1,Hr_s) -dot(s1,Hr_p);
 dot(s2,Hout_s) dot(s2,Hout_p) -dot(s2,Hr_s) -dot(s2,Hr_p) ];

Cs = [ dot(s1,Einc_s); dot(s2,Einc_s); dot(s1,Hinc_s); dot(s2,Hinc_s) ];
Cp = [ dot(s1,Einc_p); dot(s2,Einc_p); dot(s1,Hinc_p); dot(s2,Hinc_p) ];


As = F\Cs;
Ap = F\Cp ;

% Combine these 
Uin = [Einc_s, Einc_p, S_inc];
if ~tir
    P_t = [As(1)*Eout_s + As(2)*Eout_p, Ap(1)*Eout_s + Ap(2)*Eout_p, S_out] / Uin;
else
    P_t = zeros(3,3);
end
P_r = [As(3)*Er_s + As(4)*Er_p, Ap(3)*Er_s + Ap(4)*Er_p, S_r] / Uin;

% Diagnostic
Eout_for_s = As(1)*Eout_s + As(2)*Eout_p;
Hout_for_s = As(1)*Hout_s + As(2)*Hout_p;

Iin_s = dot(cross(Einc_s, conj(Hinc_s)), ada);
if ~tir
    Iout_s = dot(cross(Eout_for_s, conj(Hout_for_s)), ada);
else
    Iout_s = 0;
end

transmissionFluxRatio_s = Iout_s/Iin_s; % This should not be > 1

% MaKe sure T+R = 1
theta_i = acos(dot(k_inc,ada));
if ~tir
    theta_t = acos(dot(k_out, ada));
    T_p = (n2*cos(theta_t))/(n1*cos(theta_i))*abs(Ap(2))^2;
    T_s = (n2*cos(theta_t))/(n1*cos(theta_i))*abs(As(1))^2;
else
    T_p = 0;
    T_s = 0;
end

R_p = abs(Ap(4))^2;
R_s = abs(As(3))^2;

Test_p = T_p + R_p; 
Test_s = T_s + R_s;

[pIn, sIn] = calcPandSUnitVectors(k_inc, ada, s1);
Oin = calcO(sIn, pIn, k_inc);

%% Parallel Transport Matrix Calcs

if ~tir
    [pOut, sOut] = calcPandSUnitVectors(k_out, ada, s1);
    Oout = calcO(sOut, pOut, k_out);
    Q = Oout / Oin;
else
    pOut = zeros(3,1);
    sOut = zeros(3,1);
    Oout = eye(3);
    Q = eye(3);
end

% Section 17.2 of PLAOS
[pRef, sRef] = calcPandSUnitVectors(k_r, ada, s1);
Oref = calcO(sRef, pRef, k_r);
I_r = getIReflect;
Q_r = Oref*(I_r / Oin);
I_t = eye(3);
Q_t = Oout*(I_t / Oin);

SD_t = S_out * transpose(k_inc);
SD_r = S_r * transpose(k_inc);
transmittedActive = ~tir && ...
    isPropagatingGeometricRay(k_out, S_out, n2);
child = buildPolarizedChild(ray, hit, normal, mediumOut, ...
    Mode="isotropic", BranchType="transmitted", ...
    k=k_out, S=S_out, P=P_t, Q=Q_t, O=Oout, ...
    LocalBasis=struct('s', sOut, 'p', pOut, 'basisDirection', k_out), ...
    Index=n2, FluxNormal=ada, Active=transmittedActive, ...
    PropagationProjector=SD_t, ...
    Metadata=struct('isPropagating', transmittedActive));

children = child;

reflected = buildPolarizedChild(ray, hit, normal, mediumIn, ...
    Mode="isotropic", BranchType="reflected", ...
    k=k_r, S=S_r, P=P_r, Q=Q_r, O=Oref, ...
    LocalBasis=struct('s', sRef, 'p', pRef, 'basisDirection', k_r), ...
    Index=nr, FluxNormal=-ada, PropagationProjector=SD_r, ...
    Metadata=struct('isPropagating', true));
children = [children; reflected];

interaction.children = children;

interaction.frames.Oin = Oin;
interaction.frames.Oout = Oout;
interaction.frames.Q = Q;
interaction.frames.inputBasis = struct( ...
    's', Einc_s, ...
    'p', Einc_p, ...
    'Hs', Hinc_s, ...
    'Hp', Hinc_p, ...
    'basisDirection', k_inc);
interaction.frames.outputBasis = child.localBasis;

interaction.P = struct( ...
    'transmitted', P_t, ...
    'reflected', P_r);
interaction.Q = struct('transmitted', Q_t, ...
                       'reflected', Q_r);
interaction.coefficients = struct( ...
    'As', As, ...
    'Ap', Ap, ...
    'T_s', T_s, ...
    'T_p', T_p, ...
    'R_s', R_s, ...
    'R_p', R_p, ...
    'Test_s', Test_s, ...
    'Test_p', Test_p);
	interaction.diagnostics = struct( ...
    'F', F, ...
    'Cs', Cs, ...
    'Cp', Cp, ...
    'boundaryResidual_s', F*As - Cs, ...
    'boundaryResidual_p', F*Ap - Cp, ...
    'Iin_s', Iin_s, ...
    'Iout_s', Iout_s, ...
    'transmissionFluxRatio_s', transmissionFluxRatio_s, ...
    'Eout_s', Eout_s, ...
    'Eout_p', Eout_p, ...
    'Hout_s', Hout_s, ...
    'Hout_p', Hout_p, ...
    'tir', tir, ...
    'k_reflected', k_r, ...
    'S_reflected', S_r, ...
    'Er_s', Er_s, ...
    'Er_p', Er_p, ...
    'Hr_s', Hr_s, ...
    'Hr_p', Hr_p);

end

function n = incidentIndexFromRay(ray, mediumIn)
if isfield(ray.metadata, 'n') && ~isempty(ray.metadata.n)
    n = ray.metadata.n;
else
    n = mediumIn.IndexData{1}.n;
end
end

function tf = isPropagatingGeometricRay(k, S, n)
tf = isreal(k) && isreal(S) && isreal(n);
end



function s1 = s1FromRay(k, normal)
s1 = cross(k, normal);
if norm(s1) < 100*eps
    s1 = [1;0;0];
end
s1 = s1/norm(s1);

end
