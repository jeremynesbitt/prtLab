function interaction = traceUniaxialToUniaxial(interaction, mediumIn, mediumOut, ray, hit, normal, options)
%TRACEUNIAXIALTOISOTROPIC trace rays from uniaxial to isotropic materials
%Based on PLAOS Section 19.7.3


deg = pi/180;

k_inc = interaction.incident.k;
S_inc = interaction.incident.S;
E_inc = interaction.incident.modeE;
H_inc = interaction.incident.modeH;
iStruct = mediumIn.IndexData{1};
nO = iStruct.nO;
nE = iStruct.nE;
if strcmp(interaction.incident.mode, "ordinary")
    n_inc = iStruct.nO;
else
    n_inc = ray.metadata.n;
end
ada = interaction.normal;

Iin  = dot(cross(E_inc, conj(H_inc)), ada);
IinField = dot(cross(ray.fieldE, conj(ray.fieldH)), ada);

E_n = cross(S_inc,E_inc);
E_n = E_n/norm(E_n);

oStruct = mediumOut.IndexData{1};
n_out = oStruct.n;

% Use snell's law to calculate the exit fields.  Since this is birefringent
% to isotropic there is no issue with treating the o and e waves separately

% These two should be the same
k_out = snell_vector(k_inc, ada, n_inc, n_out);
%k_out_e = snell_vector(sout.k_tE, ada, sout.ne, 1);

%assert(norm(k_out_e-k_out_o) < sqrt(eps), 'Error!  o and e wave directions are not matched in isotropic region!');

S_out = k_out;

[E_tb, E_ta] = calcPandSUnitVectors(k_out, ada) ;
K_out = makeK(k_out);
H_ta = n_out*K_out * E_ta;
H_tb = n_out*K_out * E_tb;


%% Reflection 
% need to account for the possibility of two reflected waves for
% extraordinary wave.  

axisData = mediumIn.AxisData{1};
opticAxis = axisData.opticAxis;


a = opticAxis(:) / norm(opticAxis); % norm here is unnecessary since it is handled when sout.opticAxis is defined
% I need to construct a basis of 3 vectors, 2 orthogonal to a.  It doesn't
% matter which two vectors I choose; the a*a' accomplishes this.
epsilon = nO^2*eye(3) + (nE^2 - nO^2)*(a*a.');

% wave vector parallel to the interface.  TODO:  Update this for arbitrary ada
beta = n_inc * k_inc(1:2);

% Compute o and e reflected waves
[n_ref, k_ref, q_ref] = extraordinaryReflectedQ(beta, opticAxis, nO, nE);
qz_o = -sqrt(nO^2 - beta(1)^2 - beta(2)^2);
q_o = [beta(1); beta(2); qz_o];
k_o = q_o / nO;

K_re = makeK(k_ref);
M_re = epsilon + (n_ref*K_re)*(n_ref*K_re);
det(M_re) ;% should be close to eps
%norm((epsilon + makeK(q_ref)^2) * E_re) % should be close to eps

[Utmp,Stmp,V] = svd(M_re);
% if opticAxis, ada, and k_inc all in same direction, have a degeneracy
% that needs to be handled.
if all(opticAxis == k_inc) && all(opticAxis == ada)
   E_re = V(:,end-1);
else
   E_re = V(:,end);
end

H_re = n_ref*K_re*E_re;

K_ro = makeK(k_o);
M_ro = epsilon + (nO*K_ro)*(nO*K_ro);

det(M_ro) ; % should be close to eps

[Utmp,Stmp,V] = svd(M_ro);
%
E_ro = V(:,end);
H_ro = nO*K_ro*E_ro;

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
% This is the case III condition from chapter 19
Cn = zeros(size(Cm));


Am = F\Cm; % or F\Cs
Ap = zeros(size(Am)); % or F\Cp 

P_t = [Am(1)*E_ta+Am(2)*E_tb, zeros(size(E_ta)), S_out]*inv([E_inc, E_n, S_inc]);

% Diagnostic
E_out = Am(1)*E_ta+Am(2)*E_tb;
H_out = Am(1)*H_ta+Am(2)*H_tb;

Iout  = dot(cross(E_out, conj(H_out)), ada);
fieldE_out = P_t * ray.fieldE;
fieldScale = prtModalFieldScale(E_out, fieldE_out);
fieldH_out = fieldScale * H_out;
IoutField = dot(cross(fieldE_out, conj(fieldH_out)), ada);

Iout/Iin; % This should not be > 1

% sout_int.P_t = P_t;
% sout_int.E_tb = E_tb; % p 
% sout_int.E_ta = E_ta; % s
% sout_int.k_out = k_out; 


%% Coordinate transformation matrices (not tested)
[p1,s1] = calcPandSUnitVectors(k_inc,ada);

Oin = calcO(s1,p1,k_inc);


% o and e in same direction, so both should give same Oout
Oout   = calcO(E_ta, E_tb, k_out);

Q = Oout / Oin;

child = makeChildTemplate(ray, hit, normal, mediumOut, "transmitted");
child.k = k_out;
child.S = S_out;
child.modeE = E_out;
child.modeH = H_out;
child.fieldE = fieldE_out;
child.fieldH = fieldH_out;
child.E = child.fieldE;
child.H = child.fieldH;
child.P = P_t * ray.P;
child.Q = Q * ray.Q;
child.O = Oout;
child.localBasis = struct('s', E_ta, 'p', E_tb, 'basisDirection', k_out);
child.amplitude = norm(child.fieldE);
child.flux = real(IoutField);
child.metadata = struct( ...
    'n', n_out, ...
    'Iin', Iin, ...
    'Iout', Iout, ...
    'IinField', IinField, ...
    'IoutField', IoutField, ...
    'modalScale', fieldScale, ...
    'transmittedAmplitude', Am, ...
    'P_interface', P_t);

interaction.children = child;

interaction.frames.Oin = Oin;
interaction.frames.Oout = Oout;
interaction.frames.Q = Q;
interaction.frames.inputBasis = struct('E', E_inc, 'En', E_n, 'basisDirection', S_inc);
interaction.frames.outputBasis = child.localBasis;

interaction.P = struct('transmitted', P_t);
interaction.Q = struct('transmitted', Q);
interaction.coefficients = struct('Am', Am, 'Ap', Ap);
interaction.diagnostics = struct( ...
    'F', F, ...
    'Cm', Cm, ...
    'Cn', Cn, ...
    'boundaryResidual_m', F*Am - Cm, ...
    'boundaryResidual_n', F*Ap - Cn, ...
    'Iin', Iin, ...
    'Iout', Iout, ...
    'transmissionFluxRatio', Iout/Iin, ...
    'IinField', IinField, ...
    'IoutField', IoutField, ...
    'transmissionFluxRatioField', IoutField/IinField, ...
    'epsilon', epsilon, ...
    'opticAxis', opticAxis, ...
    'n_ref_e', n_ref, ...
    'k_ref_e', k_ref, ...
    'q_ref_e', q_ref, ...
    'k_ref_o', k_o, ...
    'E_ro', E_ro, ...
    'H_ro', H_ro, ...
    'E_re', E_re, ...
    'H_re', H_re);

end
