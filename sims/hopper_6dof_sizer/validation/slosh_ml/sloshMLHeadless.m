function out = sloshMLHeadless(R, H_tank_total, h_fill, rho, g)
% Headless replica of SLOSH-ML v1.0.0 get_geometry -> calculate ->
% compute_mechanical_analogies for a flat-bottomed right circular cylinder.
% Uses the SLOSH-ML solver functions unchanged (core, derivedparams, ...).
% Default advanced settings: N = 10, M = 5, nsteps = 200, selm = 3.
N = 10; M = 5; nsteps = 200; selm = 3;

Z = [0; H_tank_total];
r_ic = @(z) interp1(Z, [0; 0], z);
r_oc = @(z) interp1(Z, [R; R], z);
b = [0, R, 0, H_tank_total];

H_tank   = h_fill;
H_origin = h_fill;
flvol = integral(@(z) pi*(r_oc(z).^2 - r_ic(z).^2), b(3), H_origin);
z_CM  = integral(@(z) z.*(r_oc(z).^2 - r_ic(z).^2), b(3), H_origin) / ...
        integral(@(z) r_oc(z).^2 - r_ic(z).^2, b(3), H_origin);
bens  = b - [0, 0, z_CM, z_CM];
r_ic_CM = @(z) r_ic(z + z_CM);
r_oc_CM = @(z) r_oc(z + z_CM);
L   = H_tank - z_CM;
a   = r_oc_CM(L);
e   = r_ic_CM(L) / a;
L_a = L / a;
R_ic = @(Zd) r_ic_CM(Zd*a)/a;
R_oc = @(Zd) r_oc_CM(Zd*a)/a;

[c_k, lambda_k, Bzeros] = core(N, M, nsteps, bens/a, e, L_a, R_ic, R_oc, 'suppress waitbar');

intb = bens;
B = zeros(N); s_n = zeros(N,1); phi_a = zeros(N,1);
for n = 1:N
    s_n(n) = integral(@(z) (z-intb(3)/a).*R_oc(z).*phi(n,M,R_oc(z),z,Bzeros,L_a), intb(3)/a, L_a) - ...
             integral(@(z) (z-intb(3)/a).*R_ic(z).*phi(n,M,R_ic(z),z,Bzeros,L_a), intb(3)/a, L_a);
    for k = 1:N
        phi_a(n) = phi_a(n) + c_k(k,n)*phi(k,M,R_oc(L_a),L_a,Bzeros,L_a);
    end
    for m = n:N
        B(n,m) = Bint(n,m,M,Bzeros,e);
    end
end
B = triu(B) + triu(B,1).';

[~,omega_k,~,M_k,L_k,H_k,~,~,~,~,~] = derivedparams(rho,flvol,H_tank,c_k,lambda_k,g,L_a*a,a,B,s_n,phi_a,intb,selm);
[K_k, H_k_spring] = pendulum2spring(M_k, L_k, g, H_k);

out.mass      = rho*flvol;
out.omega1    = omega_k(1);        % rad/s (sorted by M_k*L_k, mode 1 first)
out.m1        = M_k(2);            % M_k(1) is the fixed mass M0
out.m0        = M_k(1);
out.L1        = L_k(2);            % pendulum length
out.H_hinge1  = H_k(2);            % pendulum hinge height (SLOSH-ML convention)
out.H_spring1 = H_k_spring(2);     % spring-mass / pendulum-bob height
out.H0        = H_k(1);            % fixed-mass height (SLOSH-ML convention)
out.K1        = K_k(2);
out.z_CM      = z_CM;
end
