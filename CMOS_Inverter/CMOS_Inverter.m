%% =========================
% CMOS Inverter with Delay and Gain
% Author: Extended Version
% =========================
clear; clc; close all;

%% =======================
% MOSFET PARAMETERS
% =======================
mu_n  = 200e-4;     % NMOS mobility (m^2/Vs)
mu_p  = 100e-4;     % PMOS mobility (m^2/Vs)
Cox   = 3e-3;       % Oxide capacitance (F/m^2)
Vth_n = 0.45;       % NMOS threshold (V)
Vth_p = -0.45;      % PMOS threshold (V)

Wn = 1e-6; Ln = 50e-9;
Wp = 2e-6; Lp = 50e-9;

beta_n = mu_n * Cox * (Wn/Ln);
beta_p = mu_p * Cox * (Wp/Lp);

Vdd = 1.2; % Supply voltage

%% ============================
% MOSFET Id functions
% ============================
Id_n = @(Vgs, Vds) ...
    (Vgs <= Vth_n) .* 0 + ...
    (Vgs > Vth_n & Vds < (Vgs - Vth_n)) .* (beta_n*((Vgs - Vth_n).*Vds - 0.5*Vds.^2)) + ...
    (Vgs > Vth_n & Vds >= (Vgs - Vth_n)) .* (0.5*beta_n*(Vgs - Vth_n).^2);

Id_p = @(Vsg, Vsd) ...
    (Vsg <= abs(Vth_p)) .* 0 + ...
    (Vsg > abs(Vth_p) & Vsd < (Vsg - abs(Vth_p))) .* (beta_p*((Vsg - abs(Vth_p)).*Vsd - 0.5*Vsd.^2)) + ...
    (Vsg > abs(Vth_p) & Vsd >= (Vsg - abs(Vth_p))) .* (0.5*beta_p*(Vsg - abs(Vth_p)).^2);

%% ============================
% Load capacitance
% ============================
CL = 10e-15; % 10 fF load capacitance

%% ============================
% CMOS Inverter VTC (Vout vs Vin)
% ============================
Vin = linspace(0, Vdd, 400);
Vout = zeros(size(Vin));

for k = 1:length(Vin)
    v = Vin(k);
    if k==1
        Vout0 = Vdd/2;
    else
        Vout0 = Vout(k-1);
    end
    
    fun = @(Vo) Id_n(v, Vo) - Id_p(Vdd - v, Vdd - Vo);
    Vout(k) = fzero(fun, Vout0);
end

figure;
plot(Vin, Vout, 'LineWidth', 2);
xlabel('V_{IN} (V)'); ylabel('V_{OUT} (V)');
title('CMOS Inverter Voltage Transfer Characteristic (VTC)');
grid on;

%% ============================
% Voltage Gain |dVout/dVin|
% ============================
gain = -diff(Vout)./diff(Vin);

figure;
plot(Vin(1:end-1), gain, 'LineWidth', 2);
xlabel('V_{IN} (V)'); ylabel('Voltage Gain |dV_{out}/dV_{in}|');
title('CMOS Inverter Voltage Gain');
grid on;

%% ============================
% Find switching point
% ============================
[Gain_max, idx_max] = max(gain);
Vsw = Vin(idx_max);   % Switching voltage
Vout_sw = Vout(idx_max);
fprintf('Switching point: Vin = %.3f V, Vout = %.3f V, Gain = %.2f\n', Vsw, Vout_sw, Gain_max);

%% ============================
% Propagation Delay Calculation (tpHL, tpLH)
% Using RC approximation: tp = CL * V/I
% ============================
% tpHL: High to Low (NMOS discharges CL)
Vout_HL = Vdd; Vin_HL = 0; % Start from Vout = Vdd
I_discharge = Id_n(Vsw, Vout_sw); % NMOS current at switching
tpHL = CL * Vdd / I_discharge;    % Propagation delay HL

% tpLH: Low to High (PMOS charges CL)
Vout_LH = 0; Vin_LH = Vdd; % Start from Vout = 0
I_charge = Id_p(Vdd - Vsw, Vout_sw); % PMOS current at switching
tpLH = CL * Vdd / I_charge;          % Propagation delay LH

fprintf('Propagation delays: tpHL = %.2e s, tpLH = %.2e s\n', tpHL, tpLH);

%% ============================
% Transient simulation (optional)
% Using simple RC charging/discharging model
% ============================
t = linspace(0, 2e-9, 1000); % 2 ns total
Vout_tr = zeros(size(t));
Vout_tr(1) = Vdd;

for k = 2:length(t)
    Vin_t = Vdd; % Step input
    I_total = Id_p(Vdd - Vin_t, Vout_tr(k-1)) - Id_n(Vin_t, Vout_tr(k-1));
    Vout_tr(k) = Vout_tr(k-1) + (I_total/CL)*(t(k)-t(k-1));
    
    % Limit Vout between 0 and Vdd
    if Vout_tr(k) > Vdd
        Vout_tr(k) = Vdd;
    elseif Vout_tr(k) < 0
        Vout_tr(k) = 0;
    end
end

figure;
plot(t*1e9, Vout_tr, 'LineWidth', 2);
xlabel('Time (ns)'); ylabel('V_{out} (V)');
title('CMOS Inverter Transient Response');
grid on;
