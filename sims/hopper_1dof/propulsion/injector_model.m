function feasible = injector_model(P_tank, Pc, dp_frac)

% Required injector ΔP
dp_required = dp_frac * Pc;

feasible = (P_tank - Pc) >= dp_required;

end
