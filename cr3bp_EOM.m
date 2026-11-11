function dydt = cr3bp_EOM(~, y, mu)
% y = [x;y;z;vx;vy;vz]
x = y(1); yy = y(2); z = y(3);
vx = y(4); vy = y(5); vz = y(6);

r1 = sqrt((x+mu)^2 + yy^2 + z^2);
r2 = sqrt((x-1+mu)^2 + yy^2 + z^2);

% partial derivatives 
Ux = x - (1-mu)*(x+mu)/r1^3 - mu*(x-1+mu)/r2^3;
Uy = yy - (1-mu)*yy/r1^3 - mu*yy/r2^3;
Uz = -(1-mu)*z/r1^3 - mu*z/r2^3;

ax = 2*vy + Ux;
ay = -2*vx + Uy;
az = Uz;

dydt = [vx; vy; vz; ax; ay; az];
end
