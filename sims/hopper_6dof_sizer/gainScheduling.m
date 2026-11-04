%% Gain Scheduling

N = 50;

% Mass profile
m_profile = linspace(OUT.Vehicle.WetMass, ...
                     OUT.Vehicle.DryMass, N);

% Preallocate gain schedule
Klqi_schedule = zeros(4, 15, N);

% Optional: store linearized system matrices
A_schedule = zeros(12, 12, N);
B_schedule = zeros(12, 4, N);

% Integral augmentation
Cpos = [eye(3), zeros(3,9)];

%% Compute gains

for i = 1:N

    % Current mass
    m = m_profile(i);

    % Hover equilibrium
    p = zeros(3,1);

    xi = [zeros(9,1); p];
    ui = [m*params.g; 0; 0; 0];

    % Linearize around current operating point
    [A,B] = HopperLinearization_lqi(xi,ui,params,m);


    Atilde = [ A,     zeros(12,3);
              -Cpos,  zeros(3,3)];

    Btilde = [B; zeros(3,4)];

    % Compute and Store gain
    Klqi_schedule(:,:,i) = lqr(Atilde, Btilde, QI, R);

    % Optional
    A_schedule(:,:,i) = A;
    B_schedule(:,:,i) = B;

end