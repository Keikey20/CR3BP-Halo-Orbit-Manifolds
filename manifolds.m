function manif = manifolds(X0, eigvals, eigvecs, mu, epsilon, numTraj, span)
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
manif = struct('type',{},'t',{},'Y',{},'X0',{});
unstable_idx = find(abs(eigvals)>1+1e-6,1,'first');
stable_idx = find(abs(eigvals)<1-1e-6,1,'first');

if isempty(unstable_idx) || isempty(stable_idx)
    warning('Could not clearly identify stable/unstable eigenvalues. Using first two eigenvectors as candidate directions.');
    unstable_idx = 1; stable_idx = 2;
end

v_un = real(eigvecs(:,unstable_idx));
v_st = real(eigvecs(:,stable_idx));
v_un = v_un / norm(v_un(1:3));
v_st = v_st / norm(v_st(1:3));

tspan_un = [0 span];
tspan_st = [0 -span]; % integrating backward in timeby negative span

for k=1:numTraj
    if mod(k,2)==1
        type = 'unstable';
        basev = v_un;
        tspan = tspan_un;
    else
        type = 'stable';
        basev = v_st;
        tspan = tspan_st;
    end
    sgn = (-1)^(k);
    amp = epsilon*(1 + 0.25*(k-1));
    Xpert = X0 + sgn*amp*basev;
    [tt, YY] = ode45(@(t,y) cr3bp_EOM(t,y,mu), linspace(tspan(1),tspan(2),ceil(abs(tspan(2)-tspan(1))/0.01)+1), Xpert, opts);
    manif(end+1).type = type; manif(end).t = tt; manif(end).Y = YY; manif(end).X0 = Xpert;
end
end
