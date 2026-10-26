function V = variantConfig(varargin)
%VARIANTCONFIG  Variant control variables for models/hopper_6dof.slx.
%   V = variantConfig()                        defaults (same as v2)
%   V = variantConfig('SLOSH', 'off', 'GROUND', 'strut', ...)
%   V = variantConfig(S)                       S.SLOSH = "off", ...
%
%   Returns VARIANT_<AXIS> = choice number for every axis in variantAxes.
%   Pass V to buildSimInput, which sets the variables on the
%   SimulationInput:
%       out = sim(buildSimInput('hopper_6dof', [], [], variantConfig('SLOSH', 'off')))
%   Unknown axes or choices are errors, so a typo can't silently run the
%   default.

A = variantAxes();
if nargin == 1 && isstruct(varargin{1})
    S = varargin{1};
else
    if mod(nargin, 2), error('variantConfig:args', 'Give axis/choice pairs, e.g. variantConfig(''SLOSH'', ''off'').'); end
    S = struct();
    for k = 1:2:nargin
        S.(upper(char(varargin{k}))) = varargin{k + 1};
    end
end
for f = fieldnames(S)'
    if ~isfield(A, f{1})
        error('variantConfig:axis', 'Unknown variant axis %s. Axes: %s.', f{1}, strjoin(fieldnames(A), ', '));
    end
end

V = struct();
for f = fieldnames(A)'
    choice = 1;
    if isfield(S, f{1})
        choice = find(A.(f{1}) == string(S.(f{1})));
        if isempty(choice)
            error('variantConfig:choice', '%s has no choice ''%s''. Choices: %s.', ...
                f{1}, string(S.(f{1})), strjoin(A.(f{1}), ', '));
        end
    end
    V.(['VARIANT_' f{1}]) = choice;
end
end
