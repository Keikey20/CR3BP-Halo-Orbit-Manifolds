function [Xcorr, Tc] = differential_correction(Xguess, Tguess, mu)
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
% initial guess
X0 = Xguess; T = Tguess;
for iter = 1:8
    [t, Y] = ode45(@(t,y) cr3bp_EOM(t,y,mu), linspace(0,T/2,400), X0, opts);
    Yh = Y(end,:)';
    fy = Yh(5); % vy at half period
    delta = 1e-7;
    Xp = X0; Xp(5) = Xp(5) + delta;
    [~, Yp] = ode45(@(t,y) cr3bp_EOM(t,y,mu), linspace(0,T/2,400), Xp, opts);
    Yhp = Yp(end,:)';
    df_dvy = (Yhp(5)-Yh(5))/delta;
    if abs(df_dvy) < 1e-12
        break
    end
    dvy = -fy/df_dvy;
    X0(5) = X0(5) + dvy;
    if abs(fy) < 1e-12
        break
    end
end
Xcorr = X0;
Tc = T;
end
