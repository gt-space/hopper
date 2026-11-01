function [GainTable] = TrajectoryFollowingGains3(Sf,Q,R,Trajectory,A,B,tgrid,unom,params,Fref)

% reduced state:
% x = [X Y Z Vx Vy Vz P Q R dth1 dth2 dth3]
% u = [T delta_p delta_y F_rcs]

Eta   = Trajectory;      % N x 12
tgrid = tgrid(:);
N     = numel(tgrid);
t0    = tgrid(1);
tf    = tgrid(end);
nx    = size(Eta,2);
nu    = size(unom,2);

%% reference state derivative from trajectory table
Etadot = zeros(size(Eta));   % N x 12

dtc = tgrid(3:end) - tgrid(1:end-2);   % (N-2)x1
Etadot(2:end-1,:) = (Eta(3:end,:) - Eta(1:end-2,:)) ./ dtc;

Etadot(1,:)   = (Eta(2,:)   - Eta(1,:))   / (tgrid(2)   - tgrid(1));
Etadot(end,:) = (Eta(end,:) - Eta(end-1,:)) / (tgrid(end) - tgrid(end-1));

%% defect: d = xref_dot - f_red(xref, uref)
d = Etadot - Fref;   % N x 12

%% solve Riccati for S backwards
[tS,Svec] = ode45(@(t,S) Sdot12(t,S,A,B,Q,R,tgrid), [tf t0], Sf(:));

tS   = flipud(tS);
Svec = flipud(Svec);

NtS = length(tS);
S_hist = zeros(nx,nx,NtS);
for k = 1:NtS
    S_hist(:,:,k) = reshape(Svec(k,:),nx,nx);
end

%% terminal condition for V
% for reduced error coordinates, terminal reference is usually zero
Vf = -Sf*(Trajectory(end,:).');

%% solve defect-based V backwards
[tV,Vvec] = ode45(@(t,V) Vdot(t,V,tS,S_hist,A,B,Q,R,d,tgrid), ...
                  [tf t0], Vf);

tV   = flipud(tV);
Vvec = flipud(Vvec);

%% interpolate S and V onto tgrid
V_on_grid = interp1(tV, Vvec, tgrid, 'linear', 'extrap').';   % nx x N

Sflat = reshape(S_hist, nx*nx, NtS).';
Sflat_on_grid = interp1(tS, Sflat, tgrid, 'linear', 'extrap');
S_on_grid = reshape(Sflat_on_grid.', nx, nx, []);

%% gains
K1 = zeros(nu,nx,N);
K2 = zeros(nu,N);

for i = 1:N
    Bi = B(:,:,i);
    Si = S_on_grid(:,:,i);
    Vi = V_on_grid(:,i);

    K1(:,:,i) = R \ (Bi' * Si);   % nu x nx
    K2(:,i)   = R \ (Bi' * Vi);   % nu x 1
end

GainTable.t  = tgrid;
GainTable.K1 = K1;
GainTable.K2 = K2;

end