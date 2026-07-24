function [n_ref, k_ref, q_ref] = extraordinaryReflectedQ(q_inc, ada, opticAxis, nO, nE)
% this func solves Fresnel's quartic equation.

% Split up k_inc into parallel and perpendicular components

qNormal = dot(q_inc,ada);
qParallel = q_inc - qNormal*ada;

a = opticAxis/norm(opticAxis);



% Now we can compute reflected waves

% qtranspose*M*q = 1 is the equation we are trying to solve

M = eye(3)/nE^2 + (1/nO^2 - 1/nE^2) * (a*a.');

A = ada.' * M * ada;
B = 2 * (ada.' * M * qParallel);
C = qParallel.' * M * qParallel - 1;

qNormalRoots = roots([A B C]);

% Need to find the correct root.  It should have energy away from the
% interface
q_ref = [];
for ii = 1:numel(qNormalRoots)
    candidate = qParallel + qNormalRoots(ii)*ada;
    energyDirection = M*candidate;

    if real(dot(energyDirection,ada)) < 0
        q_ref = candidate;
        break;
    end
end
if isempty(q_ref)
    error('extraordinaryReflectedQ:NoReflectedRoot', ...
        'Could not identify an extraordinary reflected root.');
end

n_ref = sqrt(dot(q_ref,q_ref));
k_ref = q_ref/n_ref;

end
