function [A,B] = HopperLinearization_lqi(x, u, params, mass)
% x = [X Y Z Vx Vy Vz P Q R q0 q1 q2 q3]'
% u = [T delta_p delta_y F_rcs]'   (radians)

% ----- full-state nonlinear dynamics (13x1)

h = 1e-12;                     
nx = numel(x); 
nu = numel(u);

A  = zeros(nx,nx);
B  = zeros(nx,nu);  


% ----- full 13x13 A matrix
for i = 1:nx
    dx = zeros(nx,1); 
    dx(i) = 1;
    fi = hopper_f(x + 1i*h*dx, u, params, mass);
    A(:,i) = imag(fi)/h;
end

% ----- full 13x4 B matrix
for j = 1:nu
    du = zeros(nu,1); 
    du(j) = 1;
    fj = hopper_f(x, u + 1i*h*du, params, mass);
    B(:,j) = imag(fj)/h;
end



end

function f = hopper_f(x,u,par,mass)
    X=x(1); Y=x(2); Z=x(3);
    Vx=x(4); Vy=x(5); Vz=x(6);
    P=x(7); Q=x(8); R=x(9);
    p = x(10:12);   

    T       = u(1);
    delta_p = u(2);
    delta_y = u(3);
    Tor_rcs = u(4);
    
    m = mass; 
    g = par.g; 
    I = par.I; 
    d = par.d; 

    % thrust in body
    Fx_T =  T * cos(delta_p) * cos(delta_y);
    Fy_T = -T * cos(delta_p) * sin(delta_y);
    Fz_T = -T * sin(delta_p);   
    F_Tb = [Fx_T; Fy_T; Fz_T];

    % inertial translational dynamics
    Rbn = Rbody2NED(p);
    F_n  = Rbn * F_Tb + [0;0;m*g];
    Vdot = (1/m) * F_n;
    posdot = [Vx;Vy;Vz];
    
    % rotational dynamics
    r_tvc = [-d;0;0];
    M_tvc = cross(r_tvc, F_Tb);
    M_body = [M_tvc(1)+Tor_rcs;
              M_tvc(2);
              M_tvc(3)];
    
    omega = [P;Q;R];
    omegadot = I \ (M_body - cross(omega, I*omega));
    
   
    pdot = ((1+norm(p)^2)/4)*(eye(3)+2*( (hat(p)^2+hat(p))/(1+norm(p)^2 )))*omega;
    f = [posdot;
         Vdot;
         omegadot;
         pdot];
end

function R = Rbody2NED(p)
   R=(eye(3)-hat(p))/(eye(3)+hat(p));
end


function V = hat(v)
    V = [0 -v(3) v(2);
         v(3) 0 -v(1);
         -v(2) v(1) 0];
end







