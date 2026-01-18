%% Hopper Propulsion System

function [ox_mdot, fu_mdot, tot_mdot, Pc, thrust, MR] = prop_system(thrust, ox_name, fu_name, ox_vap_P, ox_tank_P, fu_tank_P, ox_sys_CdA, fu_sys_CdA, eta_cstar, At, eps)

    %CP = py.importlib.import_module('CoolProp.CoolProp');
    
    % Hardcoded Pa - implement this as a function of altitude so it changes
    % with launch site as well
    Pa = 14.8; % psia

    % Constants
    R = 8.314; % J/mol-K
    MW = 21.683 / 1000; % g/mol --> kg/mol
    
    % Gamma set - implement CEA to calculate gamma!
    gamma = 1.26;

    % Cstar Theoretical
    cstar_theo = 1513; % m/s
    cstar_act = eta_cstar * cstar_theo; % m/s

    ox_rho = 896; % kg/m^3
    fu_rho = 786; % kg/m^2
    
    % Equations -- x(1) = ox_mdot, x(2) = fu_mdot, x(3) = Pc, x(4) = tot_mdot
    fun = @(x) [
        x(1) - (ox_sys_CdA * sqrt(2 * ox_rho * (ox_tank_P - x(3))));
        x(2) - (fu_sys_CdA * sqrt(2 * fu_rho * (fu_tank_P - x(3))));
        x(4) - (x(1) + x(2));
        cstar_act - ((x(3) * At) / x(4));
    ];
    x0 = [0.9, 0.3, 250*6894.7, 1.2]; % Initial Guess
    x = fsolve(fun, x0);
    ox_mdot  = x(1);
    fu_mdot  = x(2);
    Pc       = x(3);
    tot_mdot = x(4);

    MR = ox_mdot / fu_mdot; 
   
    
    % Solve Area-Mach Relation
    f = @(M) (eps)^2 - (1 / M^2) * (((2 / (gamma + 1)) * (1 + 0.5*(gamma - 1)*M^2))^((gamma + 1) / (gamma - 1)));
    M0 = 3;
    M = fzero(f, M0);
    
    % Pe Solve
    Pe = Pc / ((1 + 0.5*(gamma - 1)*M^2) ^ (gamma / (gamma - 1))); % Pa

    % Exit Velocity
    A = (1 / (cstar_act * ((0.5*(gamma + 1))^((-gamma - 1) / (2 * (gamma - 1))))))^2;
    T0 = gamma / (A * (R/MW)); % K 
    U_e = sqrt(((2 * gamma) / (gamma - 1)) * (R / MW) * (T0) * (1 - ((Pe / Pc)^((gamma - 1)/gamma)))); % m/s
    % import CEA tables - output all gas properties
    
    momentum_thrust = tot_mdot * U_e;
    pressure_thrust = At * eps * (Pe - Pa);
    thrust = momentum_thrust + pressure_thrust;
    % vary Pc to solve for thrust - keep constant MR

    % Assign - OX Mdot, FU mdot, Total Mdot, Pc
    

end

