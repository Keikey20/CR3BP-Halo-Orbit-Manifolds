function dydt = cr3bp_VARIATIONAL(~, Y, mu)
% Y = [state(6); Phi(36x1)]
x = Y(1); yv = Y(2); z = Y(3);
vx = Y(4); vy = Y(5); vz = Y(6);
% distances
r1 = sqrt((x+mu)^2 + yv^2 + z^2);
r2 = sqrt((x-1+mu)^2 + yv^2 + z^2);
r1_3 = r1^3; r2_3 = r2^3;
r1_5 = r1^5; r2_5 = r2^5;
% second derivatives
Uxx = 1 - (1-mu)*(1/r1_3 - 3*(x+mu)^2 / r1_5) - mu*(1/r2_3 - 3*(x-1+mu)^2 / r2_5);
Uyy = 1 - (1-mu)*(1/r1_3 - 3*(yv)^2 / r1_5) - mu*(1/r2_3 - 3*(yv)^2 / r2_5);
Uzz = -(1-mu)*(1/r1_3 - 3*(z)^2 / r1_5) - mu*(1/r2_3 - 3*(z)^2 / r2_5);
Uxy = 3*(1-mu)*(x+mu)*yv / r1_5 + 3*mu*(x-1+mu)*yv / r2_5;
Uxz = 3*(1-mu)*(x+mu)*z / r1_5 + 3*mu*(x-1+mu)*z / r2_5;
Uyz = 3*(1-mu)*yv*z / r1_5 + 3*mu*yv*z / r2_5;
% A matrix
A = zeros(6);
A(1,4)=1; A(2,5)=1; A(3,6)=1;
A(4,1)=Uxx; A(4,2)=Uxy; A(4,3)=Uxz; A(4,5)=2;
A(5,1)=Uxy; A(5,2)=Uyy; A(5,3)=Uyz; A(5,4)=-2;
A(6,1)=Uxz; A(6,2)=Uyz; A(6,3)=Uzz;
% state derivative
dstate = cr3bp_EOM(0, Y(1:6), mu);
Phi = reshape(Y(7:end),6,6);
PhiDot = A * Phi;
dydt = [dstate; reshape(PhiDot,36,1)];
end
