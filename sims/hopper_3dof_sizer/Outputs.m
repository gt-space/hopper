function OUT = Outputs(IN, VEH, TANKS, STRUCT)

%% Structures
OUT.Structures.OxTankVolume = TANKS.volumes.oxidizer;
OUT.Structures.FuelTankVolume = TANKS.volumes.fuel;

%% Vehicle Fluids
OUT.VehicleFluids.FinalMass = ...
    IN.propulsion.oxidizer_mass + IN.propulsion.fuel_mass;

OUT.VehicleFluids.MassRatio = ...
    OUT.VehicleFluids.FinalMass / ...
    (OUT.VehicleFluids.FinalMass + VEH.mass.dry);

%% Engines
OUT.Engines.Mass = VEH.mass.engine;

%% Vehicle Mass Breakdown
OUT.Vehicle.MassDistribution.Structures = STRUCT.total;
OUT.Vehicle.MassDistribution.Tanks = VEH.mass.tanks;
OUT.Vehicle.MassDistribution.Engine = VEH.mass.engine;
OUT.Vehicle.MassDistribution.Payload = IN.structures.payload_mass;

OUT.Vehicle.DryMass = VEH.mass.dry;
OUT.Vehicle.WetMass = VEH.mass.wet;

%% TWR
OUT.Vehicle.InitialTWR = ...
    IN.propulsion.nominal_thrust / ...
    (OUT.Vehicle.WetMass * IN.const.g0);

OUT.Vehicle.FinalTWR = ...
    IN.propulsion.nominal_thrust / ...
    (OUT.Vehicle.DryMass * IN.const.g0);

end
