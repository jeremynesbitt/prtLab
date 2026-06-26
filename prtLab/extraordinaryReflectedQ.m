function [n_ref, k_ref, q_ref] = extraordinaryReflectedQ(beta, opticAxis, nO, nE)
    a = opticAxis(:) / norm(opticAxis);
    bx = beta(1);
    by = beta(2);

    % q = q0 + qz*zHat
    q0 = [bx; by; 0];
    zhat = [0; 0; 1];

    M = eye(3)/nE^2 + (1/nO^2 - 1/nE^2) * (a*a.');

    A = zhat.' * M * zhat;
    B = 2 * (zhat.' * M * q0);
    C = q0.' * M * q0 - 1;

    roots_qz = roots([A B C]);

    % reflected branch: negative z phase component
    qz = roots_qz(real(roots_qz) < 0);
    qz = qz(1);

    q_ref = [bx; by; qz];
    n_ref = norm(q_ref);
    k_ref = q_ref / n_ref;
end
