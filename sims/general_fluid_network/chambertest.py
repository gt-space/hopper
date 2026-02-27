from general_fluid_network import Engine

engine = Engine(
    fuel="Kerosene",
    oxidizer="LOX",
    mdot_ox=2.6,
    mdot_fuel=1.0,
    Pc=3e6,
    eta_cstar=0.95,
    At=0.002026,
    Ae=0.01013,
    Pa=101325
)

engine.summary()