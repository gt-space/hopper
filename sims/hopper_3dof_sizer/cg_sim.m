function [cg, MoI] = cg_sim(ox_mass, fu_mass, engine_mass, engine_cg, ox_tank_mass, ox_tank_cg, ox_tank_wall_thick, fu_tank_mass, fu_tank_cg, fu_tank_wall_thick, structures_mass, structures_cg, avi_mass, avi_cg, fluids_mass, fluids_cg, ox_tank_h, fu_tank_h, tank_r, it_h)
%CG_SIM Calculates the Center of Gravity and Moment of Inertia Matrix
%   Inputs are masses (kg), CG positions (m), and dimensions (m).

    %% 1. Constants and Unit Conversions
    ox_density = 1141;   % kg/m^3
    fu_density = 800;    % kg/m^3
    al_density = 2700;   % kg/m^3 (Density of Al6061)
    
    % Engine dimensions: 4 inch radius, 1 ft length
    eng_r = 4 * 0.0254;  % inches to meters
    eng_h = 1 * 0.3048;  % ft to meters

    %% 2. End Cap Geometry & Mass (NEW)
    % Modeled as solid cylinders attached to the ends of the tank shells.
    % Height of cap = wall thickness of the respective tank.
    
    % --- Oxidizer Tank Caps ---
    h_ox_cap = ox_tank_wall_thick;
    vol_ox_cap = pi * tank_r^2 * h_ox_cap;
    m_ox_cap   = al_density * vol_ox_cap; % Mass of ONE cap
    
    % Positions (Top and Bottom of the hollow shell)
    % Z = Tank_CG +/- (Tank_H/2 + Cap_H/2)
    z_ox_cap_bot = ox_tank_cg - (ox_tank_h / 2) - (h_ox_cap / 2);
    z_ox_cap_top = ox_tank_cg + (ox_tank_h / 2) + (h_ox_cap / 2);

    % --- Fuel Tank Caps ---
    h_fu_cap = fu_tank_wall_thick;
    vol_fu_cap = pi * tank_r^2 * h_fu_cap;
    m_fu_cap   = al_density * vol_fu_cap; % Mass of ONE cap
    
    % Positions
    z_fu_cap_bot = fu_tank_cg - (fu_tank_h / 2) - (h_fu_cap / 2);
    z_fu_cap_top = fu_tank_cg + (fu_tank_h / 2) + (h_fu_cap / 2);

    %% 3. Fluid Geometry Calculations
    % Fuel Fluid
    if fu_mass > 0
        h_fu_fluid = (fu_mass / fu_density) / (pi * tank_r^2);
        % Fluid sits at the bottom of the tank shell
        z_fu_fluid = (fu_tank_cg - fu_tank_h/2) + (h_fu_fluid / 2);
    else
        h_fu_fluid = 0; z_fu_fluid = 0;
    end

    % Oxidizer Fluid
    if ox_mass > 0
        h_ox_fluid = (ox_mass / ox_density) / (pi * tank_r^2);
        z_ox_fluid = (ox_tank_cg - ox_tank_h/2) + (h_ox_fluid / 2);
    else
        h_ox_fluid = 0; z_ox_fluid = 0;
    end

    %% 4. Calculate System Center of Gravity (CG)
    % Sum of Moments = sum(mass * position)
    moment_sum = ...
        engine_mass * engine_cg + ...
        fu_tank_mass * fu_tank_cg + ...
        ox_tank_mass * ox_tank_cg + ...
        structures_mass * structures_cg + ...
        avi_mass * avi_cg + ...
        fu_mass * z_fu_fluid + ...
        ox_mass * z_ox_fluid + ...
        m_ox_cap * z_ox_cap_bot + ... 
        m_ox_cap * z_ox_cap_top + ...  
        m_fu_cap * z_fu_cap_bot + ...  
        m_fu_cap * z_fu_cap_top;      

    total_mass = engine_mass + fu_tank_mass + ox_tank_mass + ...
                 structures_mass + avi_mass + fu_mass + ox_mass + ...
                 2 * m_ox_cap + 2 * m_fu_cap;

    if total_mass > 0
        cg = moment_sum / total_mass; % return 1
    else
        cg = 0;
    end

    %% 5. Calculate Moments of Inertia (MoI)
    % Anonymous functions for inertia - shout out gemini for this one
    calc_solid_Izz = @(m, r) 0.5 * m * r^2;
    calc_solid_Ixx = @(m, r, h) (1/12) * m * (3*r^2 + h^2);
    calc_hollow_Izz = @(m, ro, ri) 0.5 * m * (ro^2 + ri^2);
    calc_hollow_Ixx = @(m, ro, ri, h) (1/12) * m * (3*(ro^2 + ri^2) + h^2);

    % --- A. Engine ---
    Izz_eng = calc_solid_Izz(engine_mass, eng_r);
    Ixx_eng = calc_solid_Ixx(engine_mass, eng_r, eng_h) + engine_mass * (engine_cg - cg)^2;

    % --- B. Fuel Tank (Hollow Shell) ---
    r_in_fu = tank_r - fu_tank_wall_thick;
    Izz_fut = calc_hollow_Izz(fu_tank_mass, tank_r, r_in_fu);
    Ixx_fut = calc_hollow_Ixx(fu_tank_mass, tank_r, r_in_fu, fu_tank_h) + fu_tank_mass * (fu_tank_cg - cg)^2;

    % --- C. Oxidizer Tank (Hollow Shell) ---
    r_in_ox = tank_r - ox_tank_wall_thick;
    Izz_oxt = calc_hollow_Izz(ox_tank_mass, tank_r, r_in_ox);
    Ixx_oxt = calc_hollow_Ixx(ox_tank_mass, tank_r, r_in_ox, ox_tank_h) + ox_tank_mass * (ox_tank_cg - cg)^2;

    % --- D. Structures --- just asssume cylinder in it
    Izz_str = calc_solid_Izz(structures_mass, tank_r);
    Ixx_str = calc_solid_Ixx(structures_mass, tank_r, it_h) + structures_mass * (structures_cg - cg)^2;

    % --- E. Avionics --- just asssume cylinder in it
    Izz_avi = calc_solid_Izz(avi_mass, tank_r);
    Ixx_avi = calc_solid_Ixx(avi_mass, tank_r, it_h) + avi_mass * (avi_cg - cg)^2;

    % --- F. Fluids ---
    if fu_mass > 0
        Izz_fuf = calc_solid_Izz(fu_mass, tank_r);
        Ixx_fuf = calc_solid_Ixx(fu_mass, tank_r, h_fu_fluid) + fu_mass * (z_fu_fluid - cg)^2;
    else
        Izz_fuf = 0; Ixx_fuf = 0;
    end

    if ox_mass > 0
        Izz_oxf = calc_solid_Izz(ox_mass, tank_r);
        Ixx_oxf = calc_solid_Ixx(ox_mass, tank_r, h_ox_fluid) + ox_mass * (z_ox_fluid - cg)^2;
    else
        Izz_oxf = 0; Ixx_oxf = 0;
    end

    % --- G. End Caps ---
    % 1. Ox Caps (Solid Cylinders)
    % Izz is same for both
    Izz_ox_cap = calc_solid_Izz(m_ox_cap, tank_r);
    Ixx_ox_cap_local = calc_solid_Ixx(m_ox_cap, tank_r, h_ox_cap);
    
    Ixx_ox_caps_total = ...
        (Ixx_ox_cap_local + m_ox_cap * (z_ox_cap_bot - cg)^2) + ... % Bottom Cap
        (Ixx_ox_cap_local + m_ox_cap * (z_ox_cap_top - cg)^2);      % Top Cap

    % 2. Fuel Caps (Solid Cylinders)
    Izz_fu_cap = calc_solid_Izz(m_fu_cap, tank_r);
    Ixx_fu_cap_local = calc_solid_Ixx(m_fu_cap, tank_r, h_fu_cap);
    
    Ixx_fu_caps_total = ...
        (Ixx_fu_cap_local + m_fu_cap * (z_fu_cap_bot - cg)^2) + ... % Bottom Cap
        (Ixx_fu_cap_local + m_fu_cap * (z_fu_cap_top - cg)^2);      % Top Cap

    %% 6. Summation
    Izz_total = Izz_eng + Izz_fut + Izz_oxt + Izz_str + Izz_avi + ...
                Izz_fuf + Izz_oxf + ...
                2*Izz_ox_cap + 2*Izz_fu_cap;
                
    Ixx_total = Ixx_eng + Ixx_fut + Ixx_oxt + Ixx_str + Ixx_avi + ...
                Ixx_fuf + Ixx_oxf + ...
                Ixx_ox_caps_total + Ixx_fu_caps_total;
    
    Iyy_total = Ixx_total; % Symmetry

    MoI = [Ixx_total, 0, 0; 
           0, Iyy_total, 0; 
           0, 0, Izz_total]; % return 2
end