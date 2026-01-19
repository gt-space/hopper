%% Hopper Propulsion System

function [ox_mdot, fu_mdot, tot_mdot, Pc, thrust, MR] = prop_system(thrust_target, ox_name, fu_name, fu_temp, ox_vap_P, ox_tank_P, fu_tank_P, ox_sys_CdA, fu_sys_CdA, eta_cstar, At, eps)
    
    %pyenv
    %CP = py.importlib.import_module('CoolProp.CoolProp');
    
    % Hardcoded Pa - implement this as a function of altitude so it changes
    % with launch site as well
    Pa = 14.8; % psia

    % Constants
    R = 8.314; % J/mol-K
    MW = 21.683 / 1000; % g/mol --> kg/mol
    
    MR_target = 3;

    % Gamma set - implement CEA to calculate gamma!
    gamma = 1.26;

    % Cstar Theoretical
    cstar_theo = 1513; % m/s
    cstar_act = eta_cstar * cstar_theo; % m/s

    ox_rho = 896; % kg/m^3
    fu_rho = 786; % kg/m^2

    % Solve Area-Mach Relation
    f = @(M) (eps)^2 - (1 / M^2) * (((2 / (gamma + 1)) * (1 + 0.5*(gamma - 1)*M^2))^((gamma + 1) / (gamma - 1)));
    M0 = 3;
    M = fzero(f, M0);
    
    % Equations -- x(1) = ox_mdot, x(2) = fu_mdot, x(3) = Pc, x(4) = tot_mdot
    fun = @(x) [
        x(7) - (sqrt(1 / ((1 / (ox_sys_CdA^2)) + (1 / (x(5)^2)))));
        x(8) - (sqrt(1 / ((1 / (fu_sys_CdA^2)) + (1 / (x(6)^2)))));
        x(1) - (x(7) * sqrt(2 * ox_rho * (ox_tank_P - x(3))));
        x(2) - (x(8) * sqrt(2 * fu_rho * (fu_tank_P - x(3))));
        x(4) - (x(1) + x(2));
        cstar_act - ((x(3) * At) / x(4));
        x(9) - (x(3) / ((1 + 0.5*(gamma - 1)*M^2) ^ (gamma / (gamma - 1))));
        x(10) - (x(1) / x(2));
        x(10) - MR_target;
        x(11) - (gamma / (((1 / (cstar_act * ((0.5*(gamma + 1))^((-gamma - 1) / (2 * (gamma - 1))))))^2) * (R/MW))); % K 
        x(12) - (sqrt(((2 * gamma) / (gamma - 1)) * (R / MW) * (x(11)) * (1 - ((x(9) / x(3))^((gamma - 1)/gamma))))); % m/s
        x(13) - (x(4) * x(12) + At * eps * (x(9) - Pa));
        x(13) - thrust_target;
    ];
    x0 = [0.9, 0.3, 250*6894.7, 1.2, 3.633E-5, 1.3144E-5, 1.8165E-5, 6.572E-6, 101325, 3, 2200, 2000, 1334]; % Initial Guess
    options = optimoptions('fsolve', 'MaxFunctionEvaluations', 8000);
    x = fsolve(fun, x0, options);
    ox_mdot  = x(1);
    fu_mdot  = x(2);
    Pc       = x(3);
    tot_mdot = x(4);
    ox_valve_CdA = x(5);
    fu_valve_CdA = x(6);
    ox_tot_CdA = x(7);
    fu_tot_CdA = x(8);
    Pe = x(9);
    MR = x(10);
    T0 = x(11);
    U_e = x(12);
    thrust = x(13);

    % import CEA tables - output all gas properties
    
    % vary Pc to solve for thrust - keep constant MR

    % Assign - OX Mdot, FU mdot, Total Mdot, Pc
    % - vary the 
    

end

