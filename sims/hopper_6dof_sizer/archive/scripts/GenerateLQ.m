clear
clc
close all

%% GENERATE TRAJECTORY

[Trajectory,unom,tgrid,m_profile,I_profile] = Trajectory3();

N = length(tgrid);

nx = 13;
nu = 4;

%% COST MATRICES

Sf = diag([1 1 1e10 1 1 1e6 1 1 1 1000 1000 1000 1000]);
Q  = diag([1 1 1e9 1 1 10 1 1 1 100 100 100 100]);
R  = diag([1000 100000 10000 10000]);

%% PARAMETERS

params = struct();

params.g = 9.81;
params.d = 1.0;
params.r_m = 0.125;

%% PREALLOCATE

A=zeros(nx,nx,N);
B=zeros(nx,nu,N);
Fref=zeros(N,nx);

%% LINEARIZATION LOOP

for i=1:N

    params.m = m_profile(i);
    params.I = I_profile(:,:,i);

    x = Trajectory(i,:)';
    u = unom(i,:)';

    [Ai,Bi,f0] = HopperLinearization_6DOF_inertial(x,u,params);
    disp("Generating " + i + "th Trim Point")
    

    A(:,:,i)=Ai;
    B(:,:,i)=Bi;
    Fref(i,:)=f0';

end

%% COMPUTE LQR
GainTable = TrajectoryFollowingGains2( ...
Sf,Q,R,Trajectory,A,B,tgrid,unom,params,Fref);

K1 = GainTable.K1;
K2grid = GainTable.K2';

%% FLATTEN K1 FOR SIMULINK

K1flat=zeros(N,52);

for i=1:N
    K1flat(i,:) = reshape(K1(:,:,i).',1,52);
    disp("Flattening " + i + "th Trim Point")
end

%% Separate figures for each K1 control input

control_labels = {'Throttle','Pitch TVC','Yaw TVC','RCS'};
state_labels   = {'X','Y','Z','Vx','Vy','Vz','P','Q','R','q0','q1','q2','q3'};

N  = length(tgrid);
nu = size(K1,1);
nx = size(K1,2);

for u = 1:nu
    figure('Name',['K1 Gains: ' control_labels{u}]);
    hold on
    for x = 1:nx
        plot(tgrid, squeeze(K1(u,x,:)), 'LineWidth', 1.5)
    end
    grid on
    xlabel('Time (s)')
    ylabel([control_labels{u} ' Gains'])
    xlim([tgrid(1) tgrid(end)])
    title(['LQR Gains for ' control_labels{u}])
    legend(state_labels, 'Location', 'eastoutside')
   % Enable data cursor mode
    dcm = datacursormode(gcf);
    set(dcm, 'Enable', 'on', 'SnapToDataVertex', 'on', ...
        'UpdateFcn', @(src,event) customCursor(event));
end

%% Custom cursor function using DisplayName
function txt = customCursor(event)
    pos = event.Position;        % [x y]
    target = event.Target;

    % Attempt to get DisplayName
    if isprop(target,'DisplayName')
        stateName = target.DisplayName;
    else
        stateName = 'Unknown';
    end

    % Build tooltip text
    txt = {['State: ' stateName], ...
           ['Time: ' num2str(pos(1),'%.3f') ' s'], ...
           ['Gain: ' num2str(pos(2),'%.6f')]};
end

disp("LQR Gains Generated")