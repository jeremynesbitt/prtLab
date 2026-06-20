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
surfaceReflective = isReflectiveSurface(mediumIn.SurfaceData{1});

% Use snell's law to calculate the exit fields.  Since this is birefringent
% to isotropic there is no issue with treating the o and e waves separately

% These two should be the same
k_out = snell_vector(k_inc, ada, n1, n2);
S_out = k_out;

[Eout_p, Eout_s] = interfaceSPBasis(k_out, ada, s1FromRay(k_inc, ada));
K_out = makeK(k_out);
Hout_s = n2 * K_out * Eout_s;
Hout_p = n2 * K_out * Eout_p;


%% Reflection 
% TODO:  Combine with 

nr = n1; 
k_r = k_inc-2*(k_inc'*ada)*ada;

[Er_p, Er_s] = interfaceSPBasis(k_r, ada, s1FromRay(k_inc, ada));

K_r = makeK(k_r);
Hr_s = nr*K_r*Er_s;
Hr_p = nr*K_r*Er_p;

% Doesn't matter whether I choose s or p here.  they point in the same
% direction
S_r = cross(Er_p,Hr_p);
S_r = S_r/norm(S_r);


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
[Einc_p, Einc_s] = interfaceSPBasis(k_inc, ada, s1);
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

if isNormalIncidence(k_inc, ada)
    t = 2*n1/(n1 + n2);
    r = (n1 - n2)/(n1 + n2);
    As = [t; 0; r; 0];
    Ap = [0; t; 0; r];
else
    As = F\Cs;
    Ap = F\Cp ;
end

% Combine these 
Uin = [Einc_s, Einc_p, S_inc];
P_t = [As(1)*Eout_s + As(2)*Eout_p, Ap(1)*Eout_s + Ap(2)*Eout_p, S_out] / Uin;
P_r = [As(3)*Er_s + As(4)*Er_p, Ap(3)*Er_s + Ap(4)*Er_p, S_r] / Uin;

% Diagnostic
Eout_for_s = As(1)*Eout_s + As(2)*Eout_p;
Hout_for_s = As(1)*Hout_s + As(2)*Hout_p;

Iin_s = dot(cross(Einc_s, conj(Hinc_s)), ada);
Iout_s = dot(cross(Eout_for_s, conj(Hout_for_s)), ada);

transmissionFluxRatio_s = Iout_s/Iin_s; % This should not be > 1

% MaKe sure T+R = 1
theta_i = acos(dot(k_inc,ada));
theta_t = acos(dot(k_out, ada));
T_p = (n2*cos(theta_t))/(n1*cos(theta_i))*abs(Ap(2))^2;
T_s = (n2*cos(theta_t))/(n1*cos(theta_i))*abs(As(1))^2;

R_p = abs(Ap(4))^2;
R_s = abs(As(3))^2;

Test_p = T_p + R_p; 
Test_s = T_s + R_s;

child = makeChildTemplate(ray, hit, normal, mediumOut, "transmitted");
child.k = k_out;
child.S = S_out;
child.fieldE = P_t * ray.fieldE;
child.fieldH = n2 * K_out * child.fieldE;
child.modeE = child.fieldE / norm(child.fieldE);
child.modeH = n2 * K_out * child.modeE;
child.E = child.fieldE;
child.H = child.fieldH;
child.P = P_t * ray.P;
[pOut, sOut] = interfaceSPBasis(k_out, ada, s1);
[pIn, sIn] = interfaceSPBasis(k_inc, ada, s1);
Oin = calcO(sIn, pIn, k_inc);
Oout = calcO(sOut, pOut, k_out);
Q = Oout / Oin;
child.Q = Q * ray.Q;
child.O = Oout;
child.localBasis = struct('s', sOut, 'p', pOut, 'basisDirection', k_out);
child.amplitude = norm(child.fieldE);
child.flux = real(dot(cross(child.fieldE, conj(child.fieldH)), ada));
child.active = ~surfaceReflective && isPropagatingGeometricRay(child.k, child.S, n2);
child.metadata = struct( ...
    'n', n2, ...
    'P_interface', P_t, ...
    'isPropagating', child.active);

children = child;

if surfaceReflective
    reflected = makeChildTemplate(ray, hit, normal, mediumIn, "reflected");
    reflected.k = k_r;
    reflected.S = S_r;
    reflected.fieldE = P_r * ray.fieldE;
    reflected.fieldH = nr * K_r * reflected.fieldE;
    reflected.modeE = reflected.fieldE / norm(reflected.fieldE);
    reflected.modeH = nr * K_r * reflected.modeE;
    reflected.E = reflected.fieldE;
    reflected.H = reflected.fieldH;
    reflected.P = P_r * ray.P;
    [pRef, sRef] = interfaceSPBasis(k_r, ada, s1);
    Oref = calcO(sRef, pRef, k_r);
    reflected.Q = (Oref / Oin) * ray.Q;
    reflected.O = Oref;
    reflected.localBasis = struct('s', sRef, 'p', pRef, 'basisDirection', k_r);
    reflected.amplitude = norm(reflected.fieldE);
    reflected.flux = abs(real(dot(cross(reflected.fieldE, conj(reflected.fieldH)), ada)));
    reflected.metadata = struct( ...
        'n', nr, ...
        'P_interface', P_r, ...
        'isPropagating', reflected.active);
    children = [children; reflected];
end

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
interaction.Q = struct('transmitted', Q);
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

function tf = isReflectiveSurface(surfaceData)
tf = false;
if isfield(surfaceData, 'isReflective')
    tf = logical(surfaceData.isReflective);
elseif isfield(surfaceData, 'isRelective')
    tf = logical(surfaceData.isRelective);
end
end

function tf = isNormalIncidence(k, normal)
tf = norm(cross(k, normal)) < 100*eps(max(1, norm(k)*norm(normal)));
end

function [p, s] = interfaceSPBasis(k, normal, fallbackS)
%INTERFACESPBASIS Build a stable p/s basis for a planar interface.

s = cross(k, normal);
if norm(s) < 100*eps
    s = fallbackS - dot(fallbackS, k)*k;
end
s = s/norm(s);
p = cross(k, s);
p = p/norm(p);

end

function s1 = s1FromRay(k, normal)
s1 = cross(k, normal);
if norm(s1) < 100*eps
    s1 = [1;0;0];
end
s1 = s1/norm(s1);

end
