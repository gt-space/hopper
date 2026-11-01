function xdot = f_inertial(x,u,g,m,l,I)
% x = [xe; ze; vx; vz; theta; q]
% u = [T; delta]
xe=x(1); ze=x(2); vx=x(3); vz=x(4); th=x(5); q=x(6);
T=u(1);  d=u(2);

phi = th + d;

xdot = zeros(6,1);
xdot(1) = vx;
xdot(2) = vz;
xdot(3) = (T/m)*cos(phi);
xdot(4) = g - (T/m)*sin(phi);
xdot(5) = q;
xdot(6) = (l/I)*T*sin(d);
end
