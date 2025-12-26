function vehicle = mass_properties(vehicle, prop)

%% Update total mass
vehicle.mass = vehicle.mass_dry + prop.m_ox + prop.m_fu;

%% Component CGs (example layout)
z_engine = vehicle.engine_z;
z_ox     = vehicle.ox_tank_z;
z_fu     = vehicle.fu_tank_z;
z_struct = vehicle.struct_cg_z;

m_engine = vehicle.engine_mass;
m_struct = vehicle.struct_mass;

%% CG
vehicle.cg_z = ...
    (m_engine*z_engine + ...
     prop.m_ox*z_ox + ...
     prop.m_fu*z_fu + ...
     m_struct*z_struct) / vehicle.mass;

%% MOI (slender-body proxy)
vehicle.Iyy = ...
    m_engine*(z_engine-vehicle.cg_z)^2 + ...
    prop.m_ox*(z_ox-vehicle.cg_z)^2 + ...
    prop.m_fu*(z_fu-vehicle.cg_z)^2 + ...
    m_struct*(z_struct-vehicle.cg_z)^2;

end
