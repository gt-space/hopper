function battery = size_batteries(IN, sim)

energy = IN.avionics.power * IN.avionics.hours;
energy = energy * IN.battery.margin;

m_batt = energy / IN.battery.energy_density;

battery.mass = m_batt;
battery.energy = energy;

end
