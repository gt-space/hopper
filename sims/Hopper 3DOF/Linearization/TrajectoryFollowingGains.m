function [GainTable] = TrajectoryFollowingGains(finalStateMaxs,stateMaxs,controlMaxs,Trajectory,A,B,tgrid)

% 3 DOF Linear quadratic ascent tracking 

%% inputs

% finalStateMaxs = 1x6 matrix of max allowable state values at terminal time
% finalStateMaxs = [xf_max,zf_max,vxf_max,vzf_max,thetaf_max,qf_max]

% stateMaxs = 1x6 matrix of max allowable state values at time t
% stateMaxs = [x_max,y_max,vx_max,vz_max,theta_max,q_max]

% controlMaxs = 1x2 matrix of max allowable control at time t
             %u1 = thrust  %u2 =tvc angle
% controlMaxs = [u1_max,u2_max];
nx = 6; nu =2;
% Trajecotry = n*5 matrix time states time history 
Eta = Trajectory;

% AX+Bu notation
% A = n*6 
% B = n*2 


% t = 1xn time history
t0 = tgrid(1); tf = tgrid(end);

%% outputs

% u1 = thrust commanded
% u2 = tvc angle commanded


%% Constraint Matricies

% 1/S_f(ii) = max acceptable value of (x_i(t_f))^2
Sf = [ (1/finalStateMaxs(1))^2 0 0 0 0 0
       0 (1/finalStateMaxs(2))^2 0 0 0 0
       0 0 (1/finalStateMaxs(3))^2 0 0 0
       0 0 0 (1/finalStateMaxs(4))^2 0 0
       0 0 0 0 (1/finalStateMaxs(5))^2 0
       0 0 0 0 0 (1/finalStateMaxs(6))^2 ];

% 1/A(ii) = (tf-t0)*max acceptable value of x_i(t)^2
Q = [  (1/stateMaxs(1))^2 0 0 0 0 0
       0 (1/stateMaxs(2))^2 0 0 0 0
       0 0 (1/stateMaxs(3))^2 0 0 0
       0 0 0 (1/stateMaxs(4))^2 0 0
       0 0 0 0 (1/stateMaxs(5))^2 0
       0 0 0 0 0 (1/stateMaxs(6))^2].*(tf-t0);

%  1/B(ii) = (tf-t0)*max acceptable value of u_i(t)^2
R = [(1/controlMaxs(1)^2) 0
     0 1/controlMaxs(2)^2].*(tf-t0);


%% Solving for S and V

%Xd is desired final state
Xd = Trajectory(end,:)';

Vf = -Sf*Xd;

[tS,Svec] = ode45(@(t,S) Sdot(t,S,A,B,Q,R,tgrid), ...
                  [tf t0], Sf(:));

%flip because we backwards integrated
Svec = flipud(Svec);
tS = flipud(tS);

NtS = length(tS);
S_hist = zeros(6,6,NtS);
for k = 1:NtS
    S_hist(:,:,k) = reshape(Svec(k,:),6,6);
end


[tV,Vvec] = ode45(@(t,V) Vdot(t,V,tS,S_hist,A,B,Q,R,Eta,tgrid), [tf t0], Vf);
Vvec = flipud(Vvec);
tV   = flipud(tV);


Sflat = reshape(S_hist, 36, NtS).';            % NtS × 36
Sflat_on_grid = interp1(tS, Sflat, tgrid);     % N × 36
S_on_grid = reshape(Sflat_on_grid.', 6, 6, []);% 6×6×N

%why transose?
V_on_grid = interp1(tV, Vvec, tgrid).';   % (N×6)' => 6×N





%% gains


% K1 = R^-1*B'*S;
% 
% K2 = R^-1*B'*V;

%why \, wrong thing inversed
N = numel(tgrid);
K1 = zeros(2,6,N);
K2 = zeros(2,N);

for i = 1:N
    Bi = B(:,:,i);              % safe: i=1..100
    Si = S_on_grid(:,:,i);
    Vi = V_on_grid(:,i);

    K1(:,:,i) = R \ (Bi' * Si);
    K2(:,i)   = R \ (Bi' * Vi);
end

GainTable.t  = tgrid;
GainTable.K1 = K1;
GainTable.K2 = K2;


end