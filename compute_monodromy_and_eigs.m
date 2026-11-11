function [PhiT, eigvals, eigvecs] = compute_monodromy_and_eigs(X0,T,mu)
Phi0 = eye(6);
Y0 = [X0; reshape(Phi0,36,1)];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
[~, Y] = ode45(@(t,y) cr3bp_VARIATIONAL(t,y,mu), linspace(0,T,2000), Y0, opts);
Yend = Y(end,:)';
PhiT = reshape(Yend(7:end),6,6);
[vecs, D] = eig(PhiT);
eigvals = diag(D);
[~, idx] = sort(abs(eigvals), 'descend');
eigvals = eigvals(idx);
eigvecs = vecs(:,idx);
end
