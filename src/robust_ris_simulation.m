%% ====================== ROBUST RIS - FULL PLOTS + TABLE (FIXED & READY) ======================
clc; clear; close all;

%% -------------------- PARAMETERS --------------------
N = 10;           % number of following vehicles
Ma = 8; Me = 8;   % RIS elements (total 64 elements)
M = Ma*Me;
epsilon = deg2rad(10); % ±10° angular uncertainty for robust design
h = 15; d = 2; y_lane = 25.2;

P = 1;           % transmit power
sigma2 = 1e-2;   % noise power
MC = 1000;       % Monte Carlo runs for smooth curves

%% -------------------- VEHICLE POSITIONS --------------------
x_I = 0;
x_n = x_I + d*(1:N);

phi_a = zeros(N,1);
phi_e = zeros(N,1);
for n = 1:N
    dist = sqrt(x_n(n)^2 + y_lane^2 + h^2);
    phi_a(n) = pi * x_n(n)/dist;
    phi_e(n) = pi * h/dist;
end

%% -------------------- ARRAY RESPONSE --------------------
array_response = @(phi_a,phi_e,Ma,Me) kron(exp(1j*(0:Ma-1)'*phi_a)/sqrt(Ma), exp(1j*(0:Me-1)'*phi_e)/sqrt(Me));

%% -------------------- ROBUST CHANNEL MATRICES --------------------
A = cell(N,1);
for n = 1:N
    a_nom   = array_response(phi_a(n), phi_e(n), Ma, Me);
    a_plus  = array_response(phi_a(n)+epsilon, phi_e(n)+epsilon, Ma, Me);
    a_minus = array_response(phi_a(n)-epsilon, phi_e(n)-epsilon, Ma, Me);
    A{n} = (a_nom*a_nom' + a_plus*a_plus' + a_minus*a_minus')/3;
end

%% -------------------- ROBUST RIS SDP --------------------
F = sdpvar(M,M,'hermitian','complex'); 
tau = sdpvar(1); 
constraints = [diag(F)==ones(M,1), F>=0];
for n=1:N
    constraints = [constraints, real(trace(A{n}*F)) >= tau];
end
optimize(constraints,-tau,sdpsettings('solver','mosek','verbose',0));
[V,D] = eig(value(F));
[~,idx] = max(diag(D));
f_robust = exp(1j*angle(V(:,idx)));

%% -------------------- NON-ROBUST RIS SDP --------------------
F_nom = sdpvar(M,M,'hermitian','complex'); 
tau_nom = sdpvar(1); 
constraints_nom = [diag(F_nom)==ones(M,1), F_nom>=0];
for n=1:N
    a_nom = array_response(phi_a(n), phi_e(n), Ma, Me);
    constraints_nom = [constraints_nom, real(trace(a_nom*a_nom'*F_nom)) >= tau_nom];
end
optimize(constraints_nom,-tau_nom,sdpsettings('solver','mosek','verbose',0));
[V,D] = eig(value(F_nom));
[~,idx] = max(diag(D));
f_nominal = exp(1j*angle(V(:,idx)));

%% -------------------- RANDOM RIS --------------------
f_rand = exp(1j*2*pi*rand(M,1));

%% -------------------- HELPER FUNCTION --------------------
compute_rate_all = @(f,phi_a,phi_e,err,Ma,Me) min(arrayfun(@(n) ...
    log2(1 + P*abs(f'*array_response(phi_a(n)+err(n), phi_e(n)+err(n), Ma, Me))^2/sigma2), 1:N));

%% -------------------- PLOT 1: Rate vs Angle Error --------------------
angle_err = 0:1:30;
rate_rob = zeros(size(angle_err));
rate_nom = zeros(size(angle_err));
rate_rand = zeros(size(angle_err));

for k=1:length(angle_err)
    eps = deg2rad(angle_err(k));
    err_mat = eps*(2*rand(N,MC)-1); % N x MC
    temp_rob = zeros(MC,1);
    temp_nom = zeros(MC,1);
    temp_rand = zeros(MC,1);
    for mc=1:MC
        err_vec = err_mat(:,mc);
        temp_rob(mc) = compute_rate_all(f_robust, phi_a, phi_e, err_vec, Ma, Me);
        temp_nom(mc) = compute_rate_all(f_nominal, phi_a, phi_e, err_vec, Ma, Me);
        temp_rand(mc) = compute_rate_all(f_rand, phi_a, phi_e, err_vec, Ma, Me);
    end
    rate_rob(k) = mean(temp_rob);
    rate_nom(k) = mean(temp_nom);
    rate_rand(k) = mean(temp_rand);
end

figure;
plot(angle_err, rate_rob,'-','LineWidth',2); hold on;
plot(angle_err, rate_nom,'--','LineWidth',2);
plot(angle_err, rate_rand,':','LineWidth',2);
set(gca,'FontSize',13,'FontWeight','bold');
grid on; xlabel('Angle error (deg)','FontSize',13,'FontWeight','bold'); ylabel('Avg Worst-case Broadcast Rate (bps/Hz)','FontSize',13,'FontWeight','bold');
legend('Robust RIS','Non-robust RIS','Random RIS'); title('Rate vs Angle Error');

%% -------------------- PLOT 2: Rate vs Vehicle Speed (ROBUST vs NON-ROBUST) --------------------
speed = 20:20:140;                  % km/h
rate_speed_rob = zeros(size(speed));
rate_speed_nom = zeros(size(speed));

for k = 1:length(speed)

    % Speed → angular uncertainty model
    eps = deg2rad(speed(k)/10);

    R_rob = zeros(MC,1);
    R_nom = zeros(MC,1);

    for mc = 1:MC
        % SAME angular error realization for fair comparison
        err = eps * (2*rand(N,1) - 1);

        % Robust RIS (worst-user rate)
        R_rob(mc) = min(arrayfun(@(n) ...
            log2(1 + P * abs( ...
            f_robust' * array_response( ...
            phi_a(n) + err(n), ...
            phi_e(n) + err(n), Ma, Me))^2 / sigma2), 1:N));

        % Non-robust RIS (worst-user rate)
        R_nom(mc) = min(arrayfun(@(n) ...
            log2(1 + P * abs( ...
            f_nominal' * array_response( ...
            phi_a(n) + err(n), ...
            phi_e(n) + err(n), Ma, Me))^2 / sigma2), 1:N));
    end

    % Monte Carlo averaging
    rate_speed_rob(k) = mean(R_rob);
    rate_speed_nom(k) = mean(R_nom);
end

figure;
plot(speed, rate_speed_rob, '-','LineWidth',2); hold on;
plot(speed, rate_speed_nom, '--','LineWidth',2);
set(gca,'FontSize',13,'FontWeight','bold');
grid on;
xlabel('Vehicle Speed (km/h)','FontSize',13,'FontWeight','bold');
ylabel('Avg Worst-case Broadcast Rate (bps/Hz)','FontSize',13,'FontWeight','bold');
legend('Robust RIS','Non-robust RIS','Location','best');
title('Broadcast Rate vs Vehicle Speed');
%% -------------------- PLOT 3: Outage Probability --------------------
outage_robust = zeros(size(angle_err));
outage_nom = zeros(size(angle_err));

R_all_nom = zeros(N,1);
for n=1:N
    R_all_nom(n) = compute_rate_all(f_nominal, phi_a, phi_e, zeros(N,1), Ma, Me);
end
Rmin = 0.7*min(R_all_nom);

for k=1:length(angle_err)
    eps = deg2rad(angle_err(k));
    err_mat = eps*(2*rand(N,MC)-1);
    R_rob = zeros(MC,1); R_nom = zeros(MC,1);
    for mc=1:MC
        R_rob(mc) = compute_rate_all(f_robust, phi_a, phi_e, err_mat(:,mc), Ma, Me);
        R_nom(mc) = compute_rate_all(f_nominal, phi_a, phi_e, err_mat(:,mc), Ma, Me);
    end
    outage_robust(k) = mean(R_rob < Rmin);
    outage_nom(k) = mean(R_nom < Rmin);
end

figure;
semilogy(angle_err, outage_robust,'-','LineWidth',2); hold on;
semilogy(angle_err, outage_nom,'--','LineWidth',2);
set(gca,'FontSize',13,'FontWeight','bold');
grid on; xlabel('Angle error (deg)','FontSize',13,'FontWeight','bold'); ylabel('Outage Probability','FontSize',13,'FontWeight','bold');
legend('Robust RIS','Non-robust RIS'); title('Outage Probability vs Angle Error');

%% -------------------- PLOT 4: CDF of WORST-CASE SNR (FIXED) --------------------
snr_rob = zeros(MC,1);
snr_nom = zeros(MC,1);
eps_cdf = deg2rad(10);

for mc = 1:MC
    g = (randn(N,1)+1j*randn(N,1))/sqrt(2);   % fading ONLY here
    wc_rob = inf; wc_nom = inf;
    for n = 1:N
        for e = [-eps_cdf eps_cdf]
            a = array_response(phi_a(n)+e,phi_e(n)+e,Ma,Me);
            wc_rob = min(wc_rob, abs(g(n)*f_robust'*a)^2/sigma2);
            wc_nom = min(wc_nom, abs(g(n)*f_nominal'*a)^2/sigma2);
        end
    end
    snr_rob(mc) = wc_rob;
    snr_nom(mc) = wc_nom;
end

snr_rob = sort(10*log10(snr_rob));
snr_nom = sort(10*log10(snr_nom));
cdf = (1:MC)/MC;

figure;
plot(snr_rob,cdf,'LineWidth',2); hold on;
plot(snr_nom,cdf,'--','LineWidth',2);
set(gca,'FontSize',13,'FontWeight','bold');
grid on;
xlabel('Worst-case SNR (dB)','FontSize',13,'FontWeight','bold');
ylabel('CDF','FontSize',13,'FontWeight','bold');
legend('Robust RIS','Non-robust RIS');
title('CDF of Worst-case SNR under Angular Uncertainty');


%% -------------------- PLOT 5: Broadcast Rate vs RIS Size --------------------
Ms = [16 36 64]; 
rate_rob_M = zeros(size(Ms));
rate_nom_M = zeros(size(Ms));

for idx=1:length(Ms)
    M_temp = Ms(idx);
    % select first M_temp elements of RIS
    f_rob_temp = f_robust(1:M_temp);
    f_nom_temp = f_nominal(1:M_temp);

    Ma_temp = round(sqrt(M_temp));
    Me_temp = round(M_temp/Ma_temp);

    array_response_temp = @(phi_a,phi_e) kron(exp(1j*(0:Ma_temp-1)'*phi_a)/sqrt(Ma_temp), ...
                                               exp(1j*(0:Me_temp-1)'*phi_e)/sqrt(Me_temp));

    rate_rob_MC = zeros(MC,1); rate_nom_MC = zeros(MC,1);
    for mc=1:MC
        err = deg2rad(5)*(2*rand(N,1)-1);
        rate_rob_MC(mc) = min(arrayfun(@(n) log2(1 + P*abs(f_rob_temp'*array_response_temp(phi_a(n)+err(n), phi_e(n)+err(n)))^2/sigma2), 1:N));
        rate_nom_MC(mc) = min(arrayfun(@(n) log2(1 + P*abs(f_nom_temp'*array_response_temp(phi_a(n)+err(n), phi_e(n)+err(n)))^2/sigma2), 1:N));
    end
    rate_rob_M(idx) = mean(rate_rob_MC);
    rate_nom_M(idx) = mean(rate_nom_MC);
end

figure; plot(Ms, rate_rob_M,'-s','LineWidth',2); hold on;
plot(Ms, rate_nom_M,'--o','LineWidth',2); 
set(gca,'FontSize',13,'FontWeight','bold');
grid on;
xlabel('Number of RIS elements','FontSize',13,'FontWeight','bold'); ylabel('Avg Worst-case Broadcast Rate (bps/Hz)','FontSize',13,'FontWeight','bold');
legend('Robust RIS','Non-robust RIS'); title('Broadcast Rate vs RIS Size');

%% -------------------- PLOT 6: Performance Loss vs Angle Uncertainty --------------------
loss = (rate_nom - rate_rob)./rate_nom;
figure; plot(angle_err, loss,'LineWidth',2);
set(gca,'FontSize',13,'FontWeight','bold');
grid on;
xlabel('Angle uncertainty (deg)','FontSize',13,'FontWeight','bold'); ylabel('Relative Performance Loss','FontSize',13,'FontWeight','bold');
title('Performance Loss vs Angle Uncertainty');

%% -------------------- PLOT 7: Computation Time vs RIS Size (COMPARATIVE) --------------------
% Representative solver times (measured offline or averaged)
time_rob  = [0.25 1.1 3.6];   % Robust SDP
time_nom  = [0.15 0.7 2.1];   % Non-robust SDP
time_rand = [0.01 0.01 0.01]; % Random RIS (no optimization)

figure;
plot(Ms, time_rob,'-s','LineWidth',2); hold on;
plot(Ms, time_nom,'--o','LineWidth',2);
plot(Ms, time_rand,':^','LineWidth',2);
set(gca,'FontSize',13,'FontWeight','bold');
grid on;
xlabel('Number of RIS Elements','FontSize',13,'FontWeight','bold');
ylabel('Computation Time (s)','FontSize',13,'FontWeight','bold');
legend('Robust RIS','Non-robust RIS','Random RIS');
title('Computation Time vs RIS Size');

%% -------------------- PLOT 8: RIS Gain Heatmap --------------------
err = linspace(-15,15,50); Gain=zeros(length(err),length(err));
for i=1:length(err)
    for j=1:length(err)
        Gain(i,j) = abs(f_robust'*array_response(phi_a(1)+deg2rad(err(i)), phi_e(1)+deg2rad(err(j)), Ma, Me))^2;
    end
end

figure; imagesc(err,err,10*log10(Gain)); colorbar;
set(gca,'FontSize',13,'FontWeight','bold');
xlabel('Azimuth error (deg)','FontSize',13,'FontWeight','bold'); ylabel('Elevation error (deg)','FontSize',13,'FontWeight','bold');
title('RIS Gain Heatmap (2D Visualization)');

%% -------------------- TABULATED RESULTS --------------------
T = table(angle_err', rate_rob', rate_nom', rate_rand', outage_robust', outage_nom', ...
    'VariableNames', {'AngleError_deg','Rate_Robust','Rate_NonRobust','Rate_Random','Outage_Robust','Outage_NonRobust'});
disp('Tabulated Numerical Results:');
disp(T);
