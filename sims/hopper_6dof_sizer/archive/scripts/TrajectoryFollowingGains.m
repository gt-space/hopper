function [GainTable] = TrajectoryFollowingGains(Sf,Q,R,Trajectory,A,B,tgrid,T_profile)

% 6 DOF Linear quadratic ascent tracking 

%x,y,z,vx,vy,vz,P,Q,R,q0,q1,q2,q3


% Trajecotry = n*5 matrix time states time history 
Eta = Trajectory;

tgrid = tgrid(:);          % N×1 no matter what comes in
N = numel(tgrid);
t0 = tgrid(1);
tf = tgrid(end);


% Eta: N×6, tgrid: N×1
Etadot = zeros(size(Eta));                % N×6

dtc = tgrid(3:end) - tgrid(1:end-2);      % (N-2)×1
Etadot(2:end-1,:) = (Eta(3:end,:) - Eta(1:end-2,:)) ./ dtc;

Etadot(1,:)   = (Eta(2,:)   - Eta(1,:))   / (tgrid(2)   - tgrid(1));
Etadot(end,:) = (Eta(end,:) - Eta(end-1,:)) / (tgrid(end) - tgrid(end-1));

% check 3dof. This was wrong
g = 9.81; m = 114; l = 0.25; I = .37;

% unom should be N×2: [T_nom, delta_nom]
% Example unom: thrust from your profile, delta = 0
unom = [T_profile(:), zeros(N,1)];   % N×2

%Fref = zeros(N,13);
%for i = 1:N
%    Fref(i,:) = f_inertial(Eta(i,:).', unom(i,:).', g, m, l, I).';
%end

%d = Etadot - Fref;   % N×6




% t = 1xn time history
t0 = tgrid(1); tf = tgrid(end);



%% Solving for S and V

%Xd is desired final state
Xd = Trajectory(end,:)';

%Vf = -Sf*Xd;

[tS,Svec] = ode45(@(t,S) Sdot(t,S,A,B,Q,R,tgrid), ...
                  [tf t0], Sf(:));

%flip because we backwards integrated
Svec = flipud(Svec);
tS = flipud(tS);

NtS = length(tS);
S_hist = zeros(13,13,NtS);
for k = 1:NtS
    S_hist(:,:,k) = reshape(Svec(k,:),13,13);
end


Vf = zeros(13,1);

%[tV,Vvec] = ode45(@(t,V) Vdot_defect(t,V,tS,S_hist,A,B,R,d,tgrid), [tf t0], Vf);

[tV,Vvec] = ode45(@(t,V) Vdot(t,V,tS,S_hist,A,B,Q,R,Eta,tgrid),[tf t0], Vf);

Vvec = flipud(Vvec);
tV   = flipud(tV);

% Interpolate V onto tgrid
V_on_grid = interp1(tV, Vvec, tgrid, 'linear', 'extrap').';   % 6×N



Sflat = reshape(S_hist, 169, NtS).';            % NtS × 169
Sflat_on_grid = interp1(tS, Sflat, tgrid);     % N × 169
S_on_grid = reshape(Sflat_on_grid.', 13, 13, []);% 13×13×N

%why transose?
%V_on_grid = interp1(tV, Vvec, tgrid).';   % (N×6)' => 6×N

%% gains

%check sizes of everything below
N = numel(tgrid);
K1 = zeros(4,13,N);
K2 = zeros(4,N);

for i = 1:N
    Bi = B(:,:,i);              % safe: i=1..100
    Si = S_on_grid(:,:,i);
    Vi = V_on_grid(:,i);

   K1(:,:,i) = (R \ (Bi' * Si));   % 2×6
K2(:,i)    = (R \ (Bi' * Vi));   % 2×1

end

GainTable.t  = tgrid;
GainTable.K1 = K1;
GainTable.K2 = K2;


end