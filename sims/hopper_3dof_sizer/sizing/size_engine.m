function engine = size_engine(IN, vehicle)

% Unpack
g0   = IN.const.g0;
Tnom = IN.propulsion.nominal_thrust;
Isp  = IN.propulsion.Isp_vac;
Pc   = IN.propulsion.chamber_pressure;
OF   = IN.propulsion.OF;

assert(~isempty(Isp) && ~isempty(Pc) && ~isempty(OF), ...
    'Engine inputs not fully defined');

%% =======================
% Mass flow
%% =======================
mdot_total = Tnom / (Isp * g0);
mdot_ox = mdot_total * OF / (1 + OF);
mdot_fu = mdot_total / (1 + OF);

%% =======================
% Chamber sizing (very rough)
%% =======================
cstar = 1500; % m/s (reasonable N2O/IPA)
At = mdot_total * cstar / Pc;

Rt = sqrt(At/pi);

%% =======================
% Structural sizing
%% =======================
sigma_allow = 2.5e8; % Pa (Al 2219)
FS = IN.const.FS_struct;

t_chamber = Pc * Rt / (sigma_allow / FS);
L_chamber = 6 * Rt;

rho = 2800; % kg/m^3
m_chamber = 2*pi*Rt*L_chamber*t_chamber*rho;

%% =======================
% Nozzle mass scaling
%% =======================
eps = 40;
m_nozzle = 0.7 * m_chamber * sqrt(eps/20);

%% =======================
% Regen penalty proxy
%% =======================
regen_factor = 1.25;

m_engine = regen_factor * (m_chamber + m_nozzle);

%% =======================
% Pack output
%% =======================
engine = struct();
engine.mdot_total = mdot_total;
engine.mdot_ox = mdot_ox;
engine.mdot_fu = mdot_fu;
engine.thrust_nom = Tnom;
engine.Isp = Isp;
engine.Pc = Pc;
engine.At = At;
engine.Rt = Rt;
engine.mass = m_engine;

end
