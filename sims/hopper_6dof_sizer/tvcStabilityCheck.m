%% TVC actuator stability check: real LQT pitch loop + actuator clamps
% Needs A_r, B_r, K1 in workspace. Avoids names A, B, t, tf (already used).
tic

kIdx = 2952;  jIn = 2;                 % schedule point, pitch gimbal input
Kk = K1(:,:,kIdx);  Ko = Kk;  Ko(jIn,:) = 0;
Lj = ss(A_r - B_r*Ko, B_r(:,jIn), Kk(jIn,:), 0);   % loop WITHOUT actuator

% {name, rate limit deg/s, accel limit deg/s^2}
acts = {'AG9113 (TVC)', 20.5, 235;
    'AG9113 @ 450', 20.5, 450;
    'AG9313 (TVC)', 13.5, 640;
    'No limits',    Inf,  Inf};
amps = [0.25 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5];   % gimbal wiggle amplitude, deg
fHz  = logspace(log10(0.05), log10(10), 60);
wv   = 2*pi*fHz;
Lr   = squeeze(freqresp(Lj, wv));  Lr = Lr(:);

nA = size(acts,1);  nAmp = numel(amps);
PM = nan(nA,nAmp); Wc = PM; GMlo = PM; GMhi = PM; DM = PM;
Nall = cell(nA,nAmp);

for a = 1:nA
    for i = 1:nAmp
        N = arrayfun(@(f) dfPoint(amps(i), f, acts{a,2}, acts{a,3}), fHz).';
        Nall{a,i} = N;
        S = allmargin(frd(Lr.*N, wv));
        if ~isempty(S.PhaseMargin)
            [PM(a,i), ix] = min(S.PhaseMargin);   Wc(a,i) = S.PMFrequency(ix);
            DM(a,i) = min(S.DelayMargin);
        end
        g = S.GainMargin;
        if any(g<1), GMlo(a,i) = 20*log10(max(g(g<1))); end
        if any(g>1), GMhi(a,i) = 20*log10(min(g(g>1))); end
    end
end

%% Table
for a = 1:nA
    fprintf('\n%s\n  Amp(deg)   PM(deg)  wc(rad/s)  GMlow(dB)  GMhigh(dB)  DelayM(ms)\n', acts{a,1});
    for i = 1:nAmp
        fprintf('  %6.2f   %7.1f   %7.2f   %8.1f   %8.1f   %8.0f\n', ...
            amps(i), PM(a,i), Wc(a,i), GMlo(a,i), GMhi(a,i), 1e3*DM(a,i));
    end
end

%% Plot 1: phase margin vs amplitude
figure; hold on; grid on
for a = 1:nA, plot(amps, PM(a,:), '-o', 'LineWidth', 1.5); end
yline(30, '--', 'PM floor 30°');
xline(0.65, ':', 'flight wiggle ~0.65°');
xlabel('Gimbal wiggle amplitude (deg)'); ylabel('Phase margin (deg)');
legend(acts(:,1), 'Location', 'southwest'); title('Phase margin vs wiggle size (real LQT loop)')

%% Plot 2: upper gain margin vs amplitude
figure; hold on; grid on
for a = 1:nA, plot(amps, GMhi(a,:), '-o', 'LineWidth', 1.5); end
yline(6, '--', 'GM floor 6 dB');
xline(0.65, ':', 'flight wiggle ~0.65°');
xlabel('Gimbal wiggle amplitude (deg)'); ylabel('Upper gain margin (dB)');
legend(acts(:,1), 'Location', 'southwest'); title('Upper gain margin vs wiggle size (real LQT loop)')

%% Plot 3: Bode overlay at one amplitude
iShow = find(amps==2);
figure
subplot(2,1,1); hold on; grid on
for a = 1:nA, semilogx(wv, 20*log10(abs(Lr.*Nall{a,iShow})), 'LineWidth', 1.5); end
yline(0,'--'); set(gca,'XScale','log'); ylabel('Magnitude (dB)')
title(sprintf('Open loop with actuator at %.1f° wiggle', amps(iShow)))
legend(acts(:,1))
subplot(2,1,2); hold on; grid on
for a = 1:nA, semilogx(wv, rad2deg(unwrap(angle(Lr.*Nall{a,iShow}))), 'LineWidth', 1.5); end
yline(-180,'--'); set(gca,'XScale','log'); ylabel('Phase (deg)'); xlabel('\omega (rad/s)')

fprintf('\nRuntime: %.1f s\n', toc);

%% Actuator describing function: servo w/ rate + accel clamps
function N = dfPoint(Amp, f, rmax, amax)
wn = 36; z = 0.7; w = 2*pi*f; T = 1/f;
dt = min(1e-3, T/2000);  tt = 0:dt:8*T-dt;
th = 0; thd = 0; y = zeros(size(tt));
for k = 1:numel(tt)
    thdd = wn^2*(Amp*sin(w*tt(k)) - th) - 2*z*wn*thd;
    thdd = max(min(thdd, amax), -amax);
    thd  = max(min(thd + thdd*dt, rmax), -rmax);
    th   = th + thd*dt;  y(k) = th;
end
m = tt >= 4*T;
s = 2*mean(y(m).*sin(w*tt(m)));  c = 2*mean(y(m).*cos(w*tt(m)));
N = complex(s, c)/Amp;
end