clear all;
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);
addpath(genpath(fullfile('..','..','_tools')));
addpath(genpath('..\..\_share\NMRtools\utils'));

%% Define the signal parameters
nt1 = 128; nt2 = 128;      % Number of points in each time dimension
nf1 = 1024; nf2 = 1024;       % Number of points in each dimension of the spectrum
% nf1 = nt1; nf2 = nt2;       % Number of points in each dimension of the spectrum
c1_ref = 500.13;                    c2_ref = 50.677748;
f1_ref = 4129.00494;                f2_ref = 5939.4320656;
c_ref = [c1_ref; c2_ref];           f_ref = [f1_ref; f2_ref];
dt1 = 3.4375e-04;                   dt2 = 6.32e-04;
fs1 = 1/dt1;                        fs2 = 1/dt2;
f1 = linspace(-fs1/2, fs1/2, nf1+1); f1 = f1(1:end-1);
f2 = linspace(-fs2/2, fs2/2, nf2+1); f2 = f2(1:end-1);
del1 = (f1 + f1_ref ) /c1_ref ;  %   Frequency scale in ppm (for even number of points)
del2 = (f2 + f2_ref ) /c2_ref ;  % 
t1 = [0:nt1-1]'*dt1;        % Time in seconds
t2 = [0:nt2-1]'*dt2;

%% Create simulated data
experiment = 2;
switch experiment
    case 0      % A single peak
        nf1 = 1024; nf2 = 1024;
        n_peaks = 1;
        true_chsh = [8.5 116];
        true_alpha = [25 10];
        true_beta = 2;
        true_intn = 1;
        true_tau = 0*[1 1]*(1e-05);
        true_theta = 0*[1 1]*(pi/3);
        t_echo = [0.0 0.0144 0.0288 0.0432 0.0576 0.072 0.0864 0.1008 0.1152 0.1296 0.1440];
        t_echo = [0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64 0.72 0.8];
    case 1      % Couples of overlapping peaks
        true_chsh = [7.02500 126.8500; 7.1000 127.0000; ...
            7.0500 116.9000; 7.1000 117.0000; ...
            7.07500 106.9500; 7.1000 107.0000; ...
            8.02500 126.8500; 8.1000 127.0000; ...
            8.0500 116.9000; 8.1000 117.0000; ...
            8.07500 106.9500; 8.1000 107.0000];
        true_alpha = [25 5; 30 3; 25 5; 30 3; 25 5; 30 3; ...
            25 5; 30 3; 25 5; 30 3; 25 5; 30 3];
        true_beta = [2 2 2 2 2 2, ...
            1 3 1 3 1 3]';
        true_intn = [0.25 0.75 0.25 0.75 0.25 0.75, ...
            0.25 0.75 0.25 0.75 0.25 0.75]';
        n_peaks = size(true_chsh, 1);    % Total number of peaks
        true_tau = 0*[1 1]*(1e-05);
        true_theta = 0*[1 1]*(pi/3);
        t_echo = [0.0144 0.0288 0.0432 0.0576 0.072 0.0864 0.1008 0.1152 0.1296 0.1440];
    case 2      % 25 random peaks
        n_peaks = 25;
        true_chsh = bsxfun(@plus, bsxfun(@times, rand(n_peaks, 2), [0.75 2.5]), [8 115]);
        true_alpha = bsxfun(@plus, bsxfun(@times, rand(n_peaks, 2), [15 10]), [20 2]);
        true_beta = 0.25 + 3*rand(n_peaks, 1);
        true_intn = 0.01 + 0.75*rand(n_peaks, 1);
        true_tau = 0*[1 1]*(1e-05);
        true_theta = 0*[1 1]*(pi/3);
        t_echo = [0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64 0.72 0.8];    %   [0.0144 0.0288 0.0432 0.0576 0.072 0.0864 0.1008 0.1152 0.1296 0.1440];    % 
end;
n_echo = numel(t_echo);

%% Construct a synthetic bicomplex signal yT
Z = getModelMatrix(true_chsh, true_alpha, [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', [0, 0]);
ZZ = transpose(conj(Z))*Z; ZZ = (ZZ + transpose(conj(ZZ)))/2;
irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
b_true = bsxfun(@times, true_intn, exp(-true_beta*t_echo));     % Inferred b
Zb = Z * b_true;
xT = reshape(Zb, nt1, nt2, n_echo);                     % Modeled signal
scale = norm(norm(xT(:,:,1)));
xT = xT / scale;    % Normalize the signal

%% Add noise
sig = 10*sqrt(8)*1.1e-04;  %        % In the 1150-1151-1155 experiments, for 8 scans, sigma is approx 1.1e-04, snr = 30.xx dB; for 1 scan, sig = sqrt(8)*1.1e-04 (1.15e-04  <->  30.5 dB; 0.003863  <->  0 dB);   x2 scans = +3 dB
nT = bicomplex(sig*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)), sig*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)));

% Compute SNR as the ratio of energies
snr = 10*log10(real(transpose(conj(xT(:)))*xT(:)) / real(transpose(conj(nT(:)))*nT(:)));
yT = xT + nT;     % Add the noise

% Compute SNR in the frequency domain for each peak separately. Take the ratio of the peak height to the std of noise
nF = real(ffthc(nT(:,:,1), [nf1, nf2], 'sinebell', 'bicomplex')) / sqrt(nt1*nt2);
sig_est = std(real(nF(:)));
zF = real(ffthc(reshape(Z/scale, nt1, nt2, []), [nf1, nf2], 'sinebell', 'bicomplex')) / sqrt(nt1*nt2);
snr_peaks = b_true(:, 1) .* squeeze(max(max(zF, [], 1), [], 2)) / sig_est;

%% Save the datasets
% save(['synth_echo_' num2str(experiment)], 'yT', 't_echo', 'c1_ref', 'c2_ref', 'f1_ref', 'f2_ref', 'dt1', 'dt2', 'true_alpha', 'true_chsh', 'true_beta', 'true_tau', 'true_theta');
% saveFID('Azara', pwd, yT, c_ref, f_ref, [dt1 dt2], 'filename', ['synth_echo_' num2str(experiment)]);
% % return;

%% Case 0.1a Study the effect of wrong parameters (only for experiment 0)
% if experiment == 0
%     maxIterRand = 3;   % 100
%     nSmpl = 25;
%     sig_arr = [sig*sqrt(64), sig*sqrt(32), sig*sqrt(16), sig*sqrt(8), sig*sqrt(4), sig*sqrt(2), sig, sig/sqrt(2), sig/sqrt(4), sig/sqrt(8), sig/sqrt(16)];
%     %sig_arr = sig/sqrt(8);
%     wrong_chsh1_arr = linspace(-0.04, 0.04, nSmpl);
%     wrong_chsh2_arr = linspace(-0.2, 0.2, nSmpl);
%     wrong_alpha1_arr = linspace(-25, 3*25, nSmpl);
%     wrong_alpha2_arr = linspace(-10, 3*10, nSmpl);
%     for j = 1 : numel(sig_arr)
%     for ii = 1 : maxIterRand
%         fprintf('Computing for sigma #%d, random case #%d of %d.\n', j, ii, maxIterRand);
%         sig = sig_arr(j);
%         nT = bicomplex(sig*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)), sig*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)));
%         snrAllTrue(j, ii) = 10*log10(real(transpose(conj(xT(:)))*xT(:)) / real(transpose(conj(nT(:)))*nT(:)));
%         yT = xT + nT;
%         yTv = reshape(yT, nt1*nt2, []);
%         for i = 1 : nSmpl
% % --------- 1. Wrong first chemical shift
%             Z = getModelMatrix(true_chsh + [wrong_chsh1_arr(i) 0], true_alpha, [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', [0, 0]);
%             irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
%             Zy = transpose(conj(Z))*yTv;
%             w = irZZ * real(Zy*conj(PHS));
%             % Fit the exponential
%             w(w < 0) = NaN;
%             indx_pos = find(w > 0);
%             linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
%             beta_hat_WP_chsh1(i, j, ii) = -linFit(1);
%             
% % --------- 2. Wrong second chemical shift
%             Z = getModelMatrix(true_chsh + [0 wrong_chsh2_arr(i)], true_alpha, [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', [0, 0]);
%             irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
%             Zy = transpose(conj(Z))*yTv;
%             w = irZZ * real(Zy*conj(PHS));
%             % Fit the exponential
%             w(w < 0) = NaN;
%             indx_pos = find(w > 0);
%             linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
%             beta_hat_WP_chsh2(i, j, ii) = -linFit(1);
%             
% % --------- 3. Wrong first alpha
%             Z = getModelMatrix(true_chsh, true_alpha + [wrong_alpha1_arr(i) 0], [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', [0, 0]);
%             irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
%             Zy = transpose(conj(Z))*yTv;
%             w = irZZ * real(Zy*conj(PHS));
%             % Fit the exponential
%             w(w < 0) = NaN;
%             indx_pos = find(w > 0);
%             linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
%             beta_hat_WP_alpha1(i, j, ii) = -linFit(1);
%             
% % --------- 4. Wrong second alpha
%             Z = getModelMatrix(true_chsh, true_alpha + [0 wrong_alpha2_arr(i)], [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', [0, 0]);
%             irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
%             Zy = transpose(conj(Z))*yTv;
%             w = irZZ * real(Zy*conj(PHS));
%             % Fit the exponential
%             w(w < 0) = NaN;
%             indx_pos = find(w > 0);
%             linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
%             beta_hat_WP_alpha2(i, j, ii) = -linFit(1);
%         end;
%     end;
%     end;
%     
% % Plot the results
% figure('Position', [340 175 1200 710]);
% subplot(2, 3, 1);
% plot(wrong_chsh1_arr*c1_ref, mean(abs(beta_hat_WP_chsh1 - true_beta), 3));
% xlim([-20 20]);
% title('^1H chemical shift', 'FontWeight', 'normal');
% xlabel('\Delta_\nu, Hz', 'FontSize', 9);
% ylim([0 0.025]);
% 
% subplot(2, 3, 2);
% plot(wrong_alpha1_arr/pi, mean(abs(beta_hat_WP_alpha1 - true_beta), 3));
% title('^1H peak width', 'FontWeight', 'normal');
% xlabel('\Delta_F_W_H_M, Hz', 'FontSize', 9);
% xlim([-8 23]);
% ylim([0 0.025]);
% 
% subplot(2, 3, 3);
% plot(mean(snrAllTrue, 2), mean(abs(beta_hat_WP_chsh1([1, 4, 7, 10, 13], :, :) - true_beta), 3)');
% title('SNR', 'FontWeight', 'normal');
% xlabel('SNR, dB', 'FontSize', 9);
% xlim([0 33]);
% ylim([0 0.15]);
% 
% subplot(2, 3, 4);
% plot(wrong_chsh2_arr*c2_ref, mean(abs(beta_hat_WP_chsh2 - true_beta), 3));
% title('^1^5N chemical shift', 'FontWeight', 'normal');
% xlabel('\Delta_\nu, Hz', 'FontSize', 9);
% xlim([-10 10]);
% ylim([0 0.025]);
% 
% subplot(2, 3, 5);
% plot(wrong_alpha2_arr/pi, mean(abs(beta_hat_WP_alpha2 - true_beta), 3));
% title('^1^5N peak width', 'FontWeight', 'normal');
% xlabel('\Delta_F_W_H_M, Hz', 'FontSize', 9);
% xlim([-3 9]);
% ylim([0 0.025]);
% 
% subplot(2, 3, 6);
% plot(mean(snrAllTrue, 2), mean(abs(beta_hat_WP_chsh2([1, 4, 7, 10, 13], :, :) - true_beta), 3)');
% title('SNR', 'FontWeight', 'normal');
% xlabel('SNR, dB', 'FontSize', 9);
% xlim([0 33]);
% ylim([0 0.15]);
%     
%     return;
% end;

%% Case 0.1b Study the effect of wrong chemical shifts wrt noise level (only for experiment 0)
if experiment == 0
    maxIterRand = 1000;
    sig_arr = 1./[4, 8, 16, 32, 64, 128, 256, 512];
    wrong_chsh1_arr = [0, true_alpha(1)/(pi*2), true_alpha(1)/pi, 1.5*true_alpha(1)/pi, 2*true_alpha(1)/pi]/c1_ref;    % linspace(-true_alpha(1)/pi, true_alpha(1)/pi, nSmpl);
    wrong_chsh2_arr = [0, true_alpha(2)/(pi*2), true_alpha(2)/pi, 1.5*true_alpha(2)/pi, 2*true_alpha(2)/pi]/c2_ref;    % linspace(-true_alpha(2)/pi, true_alpha(2)/pi, nSmpl);
    wrong_chsh1_arr = [0, 1.5*true_alpha(1)/pi]/c1_ref;    % linspace(-true_alpha(1)/pi, true_alpha(1)/pi, nSmpl);
    wrong_chsh2_arr = [0, 1.5*true_alpha(2)/pi]/c2_ref;    % linspace(-true_alpha(2)/pi, true_alpha(2)/pi, nSmpl);
    for j = 1 : numel(sig_arr)
        for ii = 1 : maxIterRand
            fprintf('Computing for sigma #%d, random case #%d of %d.\n', j, ii, maxIterRand);
            sig = sig_arr(j);
            nT = bicomplex(sig*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)), sig*(randn(nt1, nt2, n_echo) + 1i*randn(nt1, nt2, n_echo)));
            snrAllTrue(j, ii) = 10*log10(real(transpose(conj(xT(:)))*xT(:)) / real(transpose(conj(nT(:)))*nT(:)));
            yT = xT + nT;
            yTv = reshape(yT, nt1*nt2, []);
            for i = 1 : numel(wrong_chsh1_arr)
                % --------- 1. Wrong chemical shifts
                Z = getModelMatrix(true_chsh + [wrong_chsh1_arr(i) wrong_chsh2_arr(i)], true_alpha, [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', [0, 0]);
                irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
                Zy = transpose(conj(Z))*yTv;
                % Without phasing
                w = inv(transpose(conj(Z))*Z) * transpose(conj(Z))*yTv;
                w = sqrt(real(conj(w).*w));
                % Fit the exponential
                w(w < 0) = NaN;
                indx_pos = find(w > 0);
                linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
                beta_hat_WP_chsh(i, j, ii) = -linFit(1);
            end;
            
            % % % Compare with integration
%             % The peak is at (556, 474)
%             yT(:,1,:) = yT(:,1,:)/2;
%             yT(1,:,:) = yT(1,:,:)/2;
%             yF = ffthc(yT(:,:,:), [1024, 1024], 'exp', 'bicomplex') / sqrt(nt1*nt2);
%             yFr = real(yF);
%             w = reshape(sum(sum(yFr(450:650,350:600,:), 1), 2), 1, []);
%             % Fit the exponential
%             w(w < 0) = NaN;
%             indx_pos = find(w > 0);
%             linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
%             beta_hat_integr(j, ii) = -linFit(1);
            
            % % % Compare with point estimates of intensities
            yF = ffthc(yT(:,:,:), [1024, 1024], 'sinebell', 'bicomplex') / sqrt(nt1*nt2);
            yFr = real(yF);
            w = reshape(real(yF(556, 474, :)), 1, []);
            % Fit the exponential
            w(w < 0) = NaN;
            indx_pos = find(w > 0);
            linFit = polyfit(t_echo(indx_pos), log(w(indx_pos)), 1);
            beta_hat_intens(j, ii) = -linFit(1);
        end;
    end;
    
    % To compute SNR, take the ratio of the peak height to the std of noise
    nF = real(ffthc(nT(:,:,1), [nf1, nf2], 'none', 'bicomplex')) / sqrt(nt1*nt2);
    xF = real(ffthc(xT(:,:,1), [nf1, nf2], 'none', 'bicomplex')) / sqrt(nt1*nt2);
    snr_arr = max(max(xF))./sig_arr;
    % Plot the results
    figure;
    loglog(snr_arr, mean(abs(beta_hat_WP_chsh - true_beta)/true_beta, 3)');
    hold all;
    %loglog(snr_arr, mean(abs(beta_hat_integr - true_beta)/true_beta, 2)');
    loglog(snr_arr, mean(abs(beta_hat_intens - true_beta)/true_beta, 2)');
    title('SNR', 'FontWeight', 'normal');
    xlabel('SNR, dB', 'FontSize', 9);
    
    return;
end;

%% Case 0.2 Plot the results for the 25 simulated peaks experiment
if experiment == 2
    % Estimate betas
    irZZ = inv(real(ZZ)); irZZ = (irZZ + transpose(conj(irZZ)))/2;
    yTv = reshape(yT, nt1*nt2, []);
    Zy = transpose(conj(Z))*yTv;
    
    % With constrained phase for all peaks
    w = irZZ * real(Zy*conj(PHS));
    
%     % Without phasing
%     w = inv(transpose(conj(Z))*Z) * transpose(conj(Z))*yTv;
%     w = sqrt(real(conj(w).*w));
    
    w(w < 0) = NaN;
    for i = 1 : n_peaks
        indx_pos = find(w(i, :) > 0);
        linFit = polyfit(t_echo(indx_pos), log(w(i, indx_pos)), 1);     % All echoes
        beta_hat(i) = -linFit(1);
        expFit = fit(t_echo(indx_pos)',w(i, indx_pos)','exp1', 'Exclude', isnan(w(i, indx_pos)));
        beta_hat_expFit(i) = -expFit.b;
    end;
    beta_Mark = [1.95656427300000,1.08516364300000,0.525637993000000,1.81038072300000,1.70010200600000,1.63286633400000,2.51724311500000,2.11976682600000,2.87422396000000,2.09784341700000,2.84778584700000,0.769420165000000,1.08516364300000,2.47457375500000,0.615816634000000,0.474291053000000,1.17090534400000,0.253851563000000,1.50536663200000,0.291047947000000,2.32547323400000,2.98560936300000,0.525637993000000,1.63286633400000,1.30946613100000];
    beta_Mark_for_r2 = [1.95656427300000,NaN,NaN,1.81038072300000,1.70010200600000,NaN,2.51724311500000,2.11976682600000,2.87422396000000,2.09784341700000,2.84778584700000,0.769420165000000,NaN,2.47457375500000,0.615816634000000,0.474291053000000,1.17090534400000,0.253851563000000,1.50536663200000,0.291047947000000,2.32547323400000,2.98560936300000,NaN,NaN,1.30946613100000];
    beta_PINT = [1.94600000000000,0.847500000000000,0.479400000000000,1.80960000000000,1.68550000000000,3.36670000000000,2.52410000000000,2.16920000000000,2.82880000000000,2.11510000000000,2.85350000000000,0.782500000000000,1.42170000000000,2.30490000000000,0.526900000000000,0.460700000000000,1.23010000000000,0.218800000000000,1.60000000000000,0.306200000000000,2.36280000000000,3.11780000000000,1.80750000000000,0.835000000000000,1.33680000000000];
    
    % Compute r^2
    r2_true_hat = corrcoef(true_beta(~isnan(beta_Mark_for_r2)), beta_hat(~isnan(beta_Mark_for_r2)));
    r2_true_hat = r2_true_hat(1, 2)^2;
    r2_true_PINT = corrcoef(true_beta(~isnan(beta_Mark_for_r2)), beta_hat(~isnan(beta_Mark_for_r2)));
    r2_true_PINT = r2_true_PINT(1, 2)^2;
    r2_true_Mark = corrcoef(true_beta(~isnan(beta_Mark_for_r2)), beta_Mark_for_r2(~isnan(beta_Mark_for_r2)));
    r2_true_Mark = r2_true_Mark(1, 2)^2;
    r2_hat_Mark = corrcoef(beta_hat(~isnan(beta_Mark_for_r2)), beta_Mark_for_r2(~isnan(beta_Mark_for_r2)));
    r2_hat_Mark = r2_hat_Mark(1, 2)^2;
    
    
    % Compute the spectrum
    indx_plot = 1:n_peaks;    % Which peaks to mark and label
    indx_label = 1:n_peaks;
    nLines = 25;        % Number of contour lines to plot
    yF = ffthc(yT(:,:,1), [nf1, nf2], 'rcos', 'bicomplex') / sqrt(nt1*nt2);
    peaks_plot_1 = real(yF);
    peak_names = [19,23,6,1,11,10,14,12,25,20,8,17,24,18,22,2,13,16,3,4,15,21,5,9,7]';
    
    % Plot the spectrum
    figure('Position', [380 600 960 330]);
    
    subplot(1, 2, 1);
    contour(del1, del2, fliplr(rot90(peaks_plot_1, -1)), linspace(0.5*min(min(peaks_plot_1)), 0.85*max(max(peaks_plot_1)), nLines));
    hold all;
    scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 3, 'r*');
    text(true_chsh(indx_label, 1)-0.0001, true_chsh(indx_label, 2)-0.001*3.5, num2str(peak_names(indx_label(:))), 'FontSize', 7);
    set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
    ylabel('^1^5N, ppm', 'FontSize', 9); xlabel('^1H, ppm', 'FontSize', 9);
    xlim([6.5 10]); ylim([104 129]);
    xlim([7.9 8.8]);
    ylim([114.5 118]);
    
    % Plot the found - vs - true plot
    subplot(1, 2, 2);
    plot([0 4], [0 4], 'r--');
    hold all;
    %line([true_beta true_beta]', [beta_hat; beta_Mark], 'Color', 0.85*[1 1 1]);
    plot(true_beta, beta_hat, 'o', 'MarkerFaceColor', [0 0 1], 'Color', [0 0 1], 'MarkerSize', 3);
    % plot(true_beta, beta_hat_expFit, 's', 'Color', [1 0 0], 'MarkerSize', 4);
    line([true_beta true_beta]', [beta_hat; beta_hat_expFit], 'Color', 0.85*[1 1 1]);
    plot(true_beta, beta_Mark, 'd', 'MarkerFaceColor', [1 0 0], 'Color', [1 0 0], 'MarkerSize', 3);
    plot(true_beta, beta_PINT, 'v', 'MarkerFaceColor', [0 1 0], 'Color', [0 1 0], 'MarkerSize', 2.5);
    text(true_beta+0.02, beta_hat-0.02, num2str(peak_names), 'FontSize', 7);
    axis equal;
    xlim([0 3.5]); ylim([0 3.5]);
    
    %text(0.05, 3.3, ['r^2_t_r_-_m_d_l = ', num2str(r2_true_hat, '%0.3f')], 'FontSize', 9);
    %text(0.05, 3.0, ['r^2_t_r_-_i_n_g = ', num2str(r2_true_Mark, '%0.3f')], 'FontSize', 9);
    %text(0.05, 2.7, ['r^2_m_d_l_-_i_n_g = ', num2str(r2_hat_Mark, '%0.3f')], 'FontSize', 9);
    
    legend({'', 'Model-based', 'Integration'});
    xlabel('True values of R, sec^-^1', 'FontSize', 9);
    ylabel('Estimated R, sec^-^1', 'FontSize', 9);
    
    % Plot the error vs SNR.
%     snr_plot = 135.0 * true_intn / true_intn(12);      % To compute SNRs, scale according to peak intensities to make the highest intensity peak have SNR=135, as measured by Mark

    err_expFit = (beta_hat_expFit(:) - true_beta(:)) ./ true_beta(:);
    err_linFit = (beta_hat(:) - true_beta(:)) ./ true_beta(:);
    err_Mark = (beta_Mark(:) - true_beta(:)) ./ true_beta;
    err_PINT = (beta_PINT(:) - true_beta(:)) ./ true_beta;
    figure('Position', [600, 450, 750, 330]);
    plot([0, 140], [0, 0], 'r');
    hold all;
    plot(snr_peaks, err_linFit, 'bo', 'MarkerSize', 4);
    plot(snr_peaks, err_expFit, 'b.');
    plot(snr_peaks, err_Mark, 'd', 'MarkerFaceColor', [1 0 0], 'Color', [1 0 0], 'MarkerSize', 3.5);
    plot(snr_peaks, err_PINT, 'v', 'MarkerFaceColor', [0 1 0], 'Color', [0 1 0], 'MarkerSize', 3.0);
    text(snr_peaks, err_linFit, num2str(peak_names), 'FontSize', 7);
    legend({'', 'MERT - Linear fit', 'MERT - Exponential fit', 'Integration'});
    xlim([10, 115]); ylim([-1.2, 1.2]);
end;

return;

%% Plot the spectra - HC
indY = 1;     % Which plane to plot
indx_plot = 1:n_peaks;    % Which peaks to mark and label
indx_label = 1:n_peaks;
nLines = 25;        % Number of contour lines to plot
yF = ffthc(yT(:,:,indY), [nf1, nf2], 'rcos', 'bicomplex') / sqrt(nf1*nf2);
peaks_plot_1 = real(yF);
peaks_plot_2 = imag1(yF);
peaks_plot_3 = imag2(yF);
peaks_plot_4 = imag12(yF);

figure('Position', [500, 110, 870, 830]);
ax(1) = subplot(2, 2, 1);
contour(del1, del2, fliplr(rot90(peaks_plot_1, -1)), linspace(0.5*min(min(peaks_plot_1)), 0.85*max(max(peaks_plot_1)), nLines));
hold all;
scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 3, 'r*');
text(true_chsh(indx_label, 1)-0.0001, true_chsh(indx_label, 2)-0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim([6.5 10]); ylim([104 129]);
title('Real part');

ax(2) = subplot(2, 2, 2);
contour(del1, del2, fliplr(rot90(peaks_plot_2, -1)), linspace(0.5*min(min(peaks_plot_2)), 0.5*max(max(peaks_plot_2)), nLines));
hold all;
scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 3, 'r*');
text(true_chsh(indx_label, 1)-0.0001, true_chsh(indx_label, 2)-0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim([6.5 10]); ylim([104 129]);
title('Imaginary 1');

ax(3) = subplot(2, 2, 3);
contour(del1, del2, fliplr(rot90(peaks_plot_3, -1)), linspace(0.5*min(min(peaks_plot_3)), 0.5*max(max(peaks_plot_3)), nLines));
hold all;
scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 3, 'r*');
text(true_chsh(indx_label, 1)-0.0001, true_chsh(indx_label, 2)-0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim([6.5 10]); ylim([104 129]);
title('Imaginary 2');

ax(4) = subplot(2, 2, 4);
contour(del1, del2, fliplr(rot90(peaks_plot_4, -1)), linspace(0.5*min(min(peaks_plot_4)), 0.5*max(max(peaks_plot_4)), 100));
hold all;
scatter(true_chsh(indx_plot,1), true_chsh(indx_plot,2), 3, 'r*');
text(true_chsh(indx_label, 1)-0.0001, true_chsh(indx_label, 2)-0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim([6.5 10]); ylim([104 129]);
title('Imaginary 1-2');

% linkaxes(ax);
return;

%% Plot a single peak in 1D
indPlt_peak = [12];
indPlt_plane = [1 3 5 7];
% Find the indices in the plotting arrays that correspond to the chosen peak
[~, ind1] = min(abs(del1 - true_chsh(indPlt_peak,1)));
[~, ind2] = min(abs(del2 - true_chsh(indPlt_peak,2)));

yF = ffthc(yT(:,:,:), [nf1, nf2], 'none', 'bicomplex') / sqrt(nf1*nf2);
xF = ffthc(xT(:,:,:), [nf1, nf2], 'none', 'bicomplex') / sqrt(nf1*nf2);
peaks_plot_1 = real(yF);
peaks_plot_2 = real(xF);

figure('Position', [500, 100, 900, 850]);
title(['Peak #', int2str(indPlt_peak)]);
for i = 1 : 4
    pl1(i) = subplot(4,2,(i-1)*2+1);
    plot(del2, peaks_plot_1(ind1,:,indPlt_plane(i)));
    hold all
    plot(del2, peaks_plot_2(ind1,:,indPlt_plane(i)));
%     xlim(ylims);
    plot(true_chsh(indPlt_peak,2), 0, 'g*', 'MarkerSize', 1);
    ylabel(['Plane #', int2str(indPlt_plane(i))]);
    
    pl2(i) = subplot(4,2,i*2);
    plot(del1, peaks_plot_1(:,ind2,indPlt_plane(i)));
    hold all
    plot(del1, peaks_plot_2(:,ind2,indPlt_plane(i)));
%     xlim(xlims);
    plot(true_chsh(indPlt_peak,1), 0, 'g*', 'MarkerSize', 1);
end;
linkaxes(pl1); linkaxes(pl2);
return;
