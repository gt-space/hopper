function [eta] = baffleThicknessFactor(t_over_w, A_over_w)
% Baffle thickness effectiveness factor, Cole, NASA TN D-3716 (1966), Fig. 2
%
%   eta = (zeta - zeta_0)_{t/w} / (zeta - zeta_0)_{t/w = 0.04}
%
% i.e. the ratio of baffle-only damping (tank wall/tare damping removed)
% of a thick baffle to that of a thin plate, as a function of
% t/w (thickness / width) and A/w (double amplitude of fluid motion at the
% baffle edge / width). Valid for baffle Reynolds numbers 3e3 - 1.7e5.
%
% t/w <= 0.04 -> thin plate, eta = 1
% between the digitized t/w = 0.04, 0.2, 0.4 curves -> linear in t/w
% t/w > 0.4 -> held at the t/w = 0.4 curve (outside the test data)
% A/w outside each fairing's range -> held at the end value

% Fig. 2 fairings (digitized); 'nearest' extrapolation holds the end values
persistent F02 F04
if isempty(F02)
    F02 = griddedInterpolant([0.39 0.50 0.75 1.00 1.25 1.50], ...
                             [0.41 0.52 0.68 0.84 0.94 1.00], 'linear', 'nearest');
    F04 = griddedInterpolant([0.47 0.50 0.75 1.00 1.25 1.30], ...
                             [0.20 0.23 0.42 0.585 0.72 0.735], 'linear', 'nearest');
end

if t_over_w <= 0.04
    eta = ones(size(A_over_w));
    return
end

e1  = ones(size(A_over_w));
e02 = F02(A_over_w);
e04 = F04(A_over_w);

if t_over_w <= 0.04
    eta = e1;
elseif t_over_w <= 0.2
    s   = (t_over_w - 0.04) / (0.2 - 0.04);
    eta = (1 - s) * e1 + s * e02;
elseif t_over_w <= 0.4
    s   = (t_over_w - 0.2) / (0.4 - 0.2);
    eta = (1 - s) * e02 + s * e04;
else
    warning('baffleThicknessFactor:extrap', ...
        't/w = %.3f exceeds TN D-3716 data (0.4); using t/w = 0.4 curve.', t_over_w);
    eta = e04;
end

end
