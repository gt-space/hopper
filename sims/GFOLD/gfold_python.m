%% Initial Parameters
%start time
tic;

x = 0;
y = 0;
z = 55;

vx = 5;
vy = 0;
vz = -10;

xf = 50;
yf = 0;
zf = 0;

vxf = 0;
vyf = 0;
vzf = -0.1;

m0 = 113; % wetmass [kg]
Tmax = 2446; % maximum thrust [N]
specificImpulse = 225;

%% Setup for GFOLD
r0 = [x y z]; % Define initial position vector
v0 = [vx vy vz]; % Define initial velocity vector
rf = [xf yf zf]; % Define final position vector
vf = [vxf vyf vzf]; % final velocity vector
Isp = specificImpulse; % Define specific impulse

sol = gfold1(r0, v0, rf, m0, Tmax, Isp);

%% Display results
disp('Solution:');
disp(sol);

%% helper functions

function sol = gfold1(r0, v0, rf, m0, Tmax, Isp)
    pyenv('Version', 'C:\Users\joshu\AppData\Local\Programs\Python\Python311\python.exe');
    np = py.importlib.import_module('numpy');
    gfold = py.importlib.import_module('gfold');
    config = gfold.GFoldConfig();
    config.wet_mass = m0;
    config.max_thrust = Tmax;
    solver = gfold.GFoldSolver(config);

    solver.update_parameter('r_I_init', np.array(r0));
    solver.update_parameter('v_I_init', np.array(v0));
    solver.update_parameter('r_I_final', np.array(rf));
    solver.update_parameter('v_I_final', np.array(vf));

    sol = solver.solve(pyargs('verbose', true));

    
end


% end time
totalTime = toc;