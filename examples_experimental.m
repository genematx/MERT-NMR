%%
% This script presents examples of analysing 2D bicomplex NMR data in the
% time domain. Please select a specific example below.
% Required inputs are:
% 1. Time-domain relaxationdata in the Azara format
% 2. Matrices of model parameters (chemical shifts and peak widths)
%    estimated, for example during peak picking. Additionally, acquisition
%    delays in both dimensions can be supplied for propper phasing of the
%    spectrum for plotting.
% 3. An array of relaxation times t_echo.
% 
% The result of computations is the array of relation rates along with
% plots of modelled and experimntal spectra.

%% Default parameters
clear;
s = RandStream('mcg16807','Seed',0);
RandStream.setGlobalStream(s);
addpath('_tools\whittaker_smoother');   % Whittaker smoother for baseline calculation

nf1 = 1024;        nf2 = 1024;       % Number of points in each dimension of the spectrum (after zero-filling)
indY = 1;             % Use as the starting plane to plot
xlims = [6.5 10]; ylims = [104 129];    % Limits for plotting the spectra

%% Choose an experiment to run
example = 1150;
switch example
    case 1150  % T1
        load(fullfile('_data', 'fitted_pars_1150_108_ph'));
        fid = fopen(fullfile('_data', '1150my_t1t2.spc'));
        t_echo = [0.08 0.16 0.24 0.4 0.56 0.64 0.72 0.8];
    case 1151  % T1rho
        load(fullfile('_data', 'fitted_pars_1151_108_ph'));
        fid = fopen(fullfile('_data', '1151my_t1t2.spc'));
        t_echo = [0.001 0.005 0.01 0.018 0.035 0.055 0.08 0.11 0.15 0.2];
    case 1155  % T2
        load(fullfile('_data', 'fitted_pars_1155_108_ph'));
        fid = fopen(fullfile('_data', '1155my_t1t2.spc'));
        t_echo = [0.0144 0.0288 0.0432 0.0576 0.072 0.0864 0.1008 0.1152 0.1296 0.144];
end;

% % Uncomment to use the original parameters determined by peak-picking without optimization
% best_chsh(1:size(orig_chsh,1), :) = orig_chsh;
% best_alpha(1:size(orig_alpha, 1), :) = orig_alpha;

% % Only analyze the first 104 peaks
% best_chsh = best_chsh(1:104, :);
% best_alpha = best_alpha(1:104, :);

%% Load the data, define parameters of the spectrum (spectrometer frequencies, offsets, dwell times)
%% Load the data
procdat = fread(fid,inf,'float',0,'l');
fclose(fid);
dim1 = 1024; %direct dimension points (we get these for "free")
dim2 = 256;  %indirect dimension points (these we would want to undersample)
dim3 = numel(t_echo);    %time delays for relaxation fit
block1 = 128;
block2 = 16;
block3 = 2;

nb1 = ceil(dim1/block1);%floor(length(procdat)/dim2/block1);
nb2 = ceil(dim2/block2);%floor(length(procdat)/(nb1*block1))/block2);
nb3 = ceil(dim3/block3);%floor(length(procdat)/(nb1*block1))/block2);
temp = reshape(procdat,[block1,block2,block3,nb1,nb2,nb3]);
temp = permute(temp,[1 4 2 5 3 6]);
temp = reshape(temp,[nb1*block1,nb2*block2,nb3*block3]);

yT = bicomplex((temp(1:2:end, 2:2:end, :)+1i*temp(2:2:end, 2:2:end, :)), ...
    -(temp(1:2:end, 1:2:end, :)+1i*temp(2:2:end, 1:2:end, :)));
yT = yT(1:440, :, :);               % Cut the part left due to the initial delay
yT = yT / norm(norm(yT(:,:,1)));
% Phase for plotting (only for the T1rho and T2 experiments)
if example == 1151
    yT = yT .* bicomplex( exp(1i*0)*ones(size(yT)), exp(1i*pi)*ones(size(yT)), ones(size(yT)) );
end;
if example == 1155
    yT = yT .* bicomplex( exp(1i*0)*ones(size(yT)), exp(1i*pi/2)*ones(size(yT)), ones(size(yT)) );
end;
%% Define the signal parameters
nt = size(yT); nt1 = nt(1); nt2 = nt(2);      % Number of points in each time dimension
% nf1 = nt1; nf2 = nt2;       % Number of points in each dimension of the spectrum
c1_ref = 500.13;                    c2_ref = 50.677748;
f1_ref = 2360;                      f2_ref = 6081.32976;
f1_ref = f1_ref + 0.038*c1_ref;     f2_ref = f2_ref - 2.8*c2_ref;     % Correction to make the scales match with the provided analyzed data
dt1 = 1e-04;                        dt2 = 6.32e-04;
fs1 = 1/dt1;                        fs2 = 1/dt2;
f1 = linspace(-fs1/2, fs1/2, nf1+1); f1 = f1(1:end-1);
f2 = linspace(-fs2/2, fs2/2, nf2+1); f2 = f2(1:end-1);
del1 = (f1 + f1_ref ) /c1_ref ;  %   Frequency scale in ppm (for even number of points)
del2 = (f2 + f2_ref ) /c2_ref ;  % 
t1 = [0:nt1-1]'*dt1;           % Time in seconds
t2 = [0:nt2-1]'*dt2;
n_echo = numel(t_echo);
%% Cut only an interesting part of the spectrum (i.e. 6.2...10.2 ppm in the 1H direction)
ind_cut = [234:361];      % Indices in the frequency domain to cut (an even number)          ind_cut = [80:207];
yt1 = NaN(nt1, 2*nt2, n_echo);     % Data in time domain, complex values along the 1st dimension, R-I-R-I-... along the second
yt1(:, 1:2:end, :) = cmpl1(yT(:,:,:)); yt1(:, 2:2:end, :) = cmpl2(yT(:,:,:));
yf1t2 = fftshift(fft(yt1, nt1, 1), 1);
yf1t2 = yf1t2(ind_cut, :, :);
f1 = linspace(-fs1/2, fs1/2, nt1+1); f1 = f1(1:end-1);  % Frequency vector for the computed FT
yt1t2 = ifft(ifftshift(yf1t2, 1), [], 1);   % Back to the time domain
yT = bicomplex(yt1t2(:, 1:2:end, :), yt1t2(:, 2:2:end, :));

% Redefine the parameters for the first dimension
nt1 = numel(ind_cut);      % New number of time points
c1_ref = c1_ref;   % Stays the same
f1 = f1(ind_cut) + f1_ref;
f1_ref = f1(nt1/2+1);   % Set the middle point in the spectrum to be the reference frequency
f1 = f1 - f1_ref;  % New frequency range (in Hz)
fs1 = 2*abs(f1(1));     % Sampling frequency
dt1 = 1/fs1;       % Dwell time
% Generate new range vector
f1 = linspace(-fs1/2, fs1/2, nf1+1); f1 = f1(1:end-1);
del1 = (f1 + f1_ref) /c1_ref ;  %   Frequency scale in ppm (for even number of points)

t1 = [0:nt1-1]'*dt1;        % Time in seconds
[tx, ty] = meshgrid(t2, t1);
[fx, fy] = meshgrid(del2, del1);

%% Estimate noise level
yFr = real(ffthc(yT(:,:,indY), [nt1, nt2], 'none', 'bicomplex')) / sqrt(nt1*nt2);
f1n = linspace(-fs1/2, fs1/2, nt1+1); f1n = f1n(1:end-1);
f2n = linspace(-fs2/2, fs2/2, nt2+1); f2n = f2n(1:end-1);
del1n = (f1n + f1_ref ) /c1_ref ;  %   Frequency scale in ppm (for even number of points)
del2n = (f2n + f2_ref ) /c2_ref ;
for i = 1 : size(best_chsh, 1)
    yFr(and(del1n > best_chsh(i, 1)-0.5, del1n < best_chsh(i, 1)+0.5), :) = NaN;
    yFr(:, and(del2n > best_chsh(i, 2)-0.5, del2n < best_chsh(i, 2)+0.5)) = NaN;
end
yFr_bsln = wsmooth(yFr,linspace(0,1,nt2),linspace(0,1,nt1), 2.5);     % Aproximate the baseline
nFr = yFr - yFr_bsln;              % Only noise
nFr = nFr(:);
nFr = nFr(~isnan(nFr));
sig_est = std(nFr);

%% Fit the baseline
% Construct the model signal
Z = getModelMatrix(best_chsh, best_alpha, [c1_ref c2_ref]', [f1_ref f2_ref]', t1, t2, 'tau', best_tau);
ZZ = transpose(conj(Z))*Z;
Zy = transpose(conj(Z)) * reshape(yT(:, :, indY), [], 1);
ZZi = inv(ZZ);
b = ZZi*Zy;     % Inferred vector of intensities b

xT = reshape(Z * b, nt1, nt2);                     % Ideal modeled signal in the time domain
yF = ffthc(yT(:,:,indY), [nt1, nt2], 'none', 'bicomplex');     % Measured spectrum
xF = ffthc(xT, [nt1, nt2], 'none', 'bicomplex');               % Modelled spectrum
rF = yF - xF;                                                  % Residual spectrum

% Compute the baseline by smoothing the residual
% % Case 1. Smoothen in both dimensions
% bF = bicomplex((wsmooth(real(rF)) + 1i*wsmooth(imag1(rF))), ...
%     (wsmooth(imag2(rF)) + 1i*wsmooth(imag12(rF))) );
% Case 2. Smoothen only in the direct dimension
bF = bicomplex(zeros(nt1, nt2), zeros(nt1, nt2) );
for i = 1 : nt2
    bF(:, i) = bicomplex((wsmooth(real(rF(:, i))) + 1i*wsmooth(imag1(rF(:, i)))), ...
    (wsmooth(imag2(rF(:, i))) + 1i*wsmooth(imag12(rF(:, i)))) );
end;
bT = iffthc(bF, [nt1, nt2]);
% % Case 3. Constant baseline
% bT = bicomplex([1; zeros(nt1*nt2-1, 1)], zeros(nt1*nt2, 1));

% Add the baseline to the model matrix, renormalize and update ZZ
Z = [Z bT(:)];
Z = Z*diag(1./sqrt(diag(real(transpose(conj(Z))*Z))));     % Normalize Z
ZZ = transpose(conj(Z))*Z;
ZZi = inv(ZZ);

%% Construct the model signal
ZY = transpose(conj(Z)) * reshape(yT, [], n_echo);
b = ZZi*ZY;     % Inferred b
Zb = Z * b;                  % Plot WITH baseline
%    Zb = Z(:, 1:n_peaks) * b(1:n_peaks,:);             % Plot WITHOUT baseline
xT = reshape(Zb, nt1, nt2, n_echo);                     % Modeled signal

%% Compute spectra for plotting
% % Absolute values of bicomplex data
%peaks_plot_1 = abs(modc(ffthc(yT(:,:,:), [nf1, nf2], 'hann', 'bicomplex')));
%peaks_plot_2 = abs(modc(ffthc(xT(:,:,:), [nf1, nf2], 'hann', 'bicomplex')));

% Real spectra of reduced complex data
peaks_plot_1 = real(ffthc(yT(:,:,:), [nf1, nf2], 'sinebell', 'complex'));
peaks_plot_2 = real(ffthc(xT(:,:,:), [nf1, nf2], 'sinebell', 'complex'));
%% Plot the spectra
n_peaks = size(best_chsh, 1);    % Total number of peaks
indx_plot = 1:n_peaks;    % Which peaks to mark and label
indx_label = 1:n_peaks;

lvl_min = 0.04*max([peaks_plot_1(:); peaks_plot_2(:)]);
lvl_max = 0.85*max([peaks_plot_1(:); peaks_plot_2(:)]);

figure('Position', [180, 120, 1050, 530]);
ax(1) = subplot(1, 2, 1);
grid on; hold all;
contour(del1, del2, fliplr(rot90(peaks_plot_1(:,:,indY), -1)), linspace(lvl_min, lvl_max, 8));
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 3, 'r*');
text(best_chsh(indx_label, 1)+0.001*10, best_chsh(indx_label, 2)+0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim(xlims); ylim(ylims);

title('Observed spectrum');

ax(2) = subplot(1, 2, 2);
grid on; hold all;
contour(del1, del2, fliplr(rot90(peaks_plot_2(:,:,indY), -1)), linspace(lvl_min, lvl_max, 8));
line([best_chsh(indx_plot,1) best_chsh(indx_plot,1)]', [best_chsh(indx_plot,2) best_chsh(indx_plot,2)]');
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 5, 'ro');
scatter(best_chsh(indx_plot,1), best_chsh(indx_plot,2), 3, 'r*');
text(best_chsh(indx_label, 1)+0.001*10, best_chsh(indx_label, 2)+0.001*3.5, num2str(indx_label(:)), 'FontSize', 7);
set(gca,'XDir','reverse'); set(gca,'YDir','reverse');
ylabel('^1^5N, ppm'); xlabel('^1H, ppm');
xlim(xlims); ylim(ylims);
title('Modeled spectrum');

linkaxes(ax);
% return;

%% Plot a single peak in 1D
indPlt_peak = 37;            % Choose which peak to plot
indPlt_plane = [1 3 5 7];    % Choose which planes to show
% Find the indices in the plotting arrays that correspond to the chosen peak
[~, ind1] = min(abs(del1 - best_chsh(indPlt_peak,1)));
[~, ind2] = min(abs(del2 - best_chsh(indPlt_peak,2)));

figure('Position', [500, 100, 900, 850]);
title(['Peak #', int2str(indPlt_peak)]);
for i = 1 : 4
    pl1(i) = subplot(4,2,(i-1)*2+1);
    plot(del2, peaks_plot_1(ind1,:,indPlt_plane(i)));
    hold all
    plot(del2, peaks_plot_2(ind1,:,indPlt_plane(i)));
    xlim(ylims);
    ylim([-5, 65]);
    plot(best_chsh(indPlt_peak,2), 0, 'g*', 'MarkerSize', 1);
    ylabel(['Plane #', int2str(indPlt_plane(i))]);
    if i == 1
        legend('Measured data', 'Fitted model')
    end;
    
    pl2(i) = subplot(4,2,i*2);
    plot(del1, peaks_plot_1(:,ind2,indPlt_plane(i)));
    hold all
    plot(del1, peaks_plot_2(:,ind2,indPlt_plane(i)));
    xlim(xlims);
    ylim([-5, 65]);
    plot(best_chsh(indPlt_peak,1), 0, 'g*', 'MarkerSize', 1);
end;
linkaxes(pl1); linkaxes(pl2);

% % Plot residuals
% figure('Position', [500, 100, 900, 850]);
% title(['Peak #', int2str(indPlt_peak)]);
% for i = 1 : 4
%     pl1(i) = subplot(4,2,(i-1)*2+1);
%     plot(del2, peaks_plot_1(ind1,:,indPlt_plane(i)) - peaks_plot_2(ind1,:,indPlt_plane(i)));
%     hold all
%     xlim(ylims);
%     ylim([-4.5, 4.5]);
%     plot(best_chsh(indPlt_peak,2), 0, 'g*', 'MarkerSize', 1);
%     ylabel(['Plane #', int2str(indPlt_plane(i))]);
%     
%     pl2(i) = subplot(4,2,i*2);
%     plot(del1, peaks_plot_1(:,ind2,indPlt_plane(i)) - peaks_plot_2(:,ind2,indPlt_plane(i)));
%     hold all
%     xlim(xlims);
%     ylim([-4.5, 4.5]);
%     plot(best_chsh(indPlt_peak,1), 0, 'g*', 'MarkerSize', 1);
% end;
% linkaxes(pl1); linkaxes(pl2);
% %return;

%% Estimate T1/T2/T1rho
% Estimate peak intensities
w = ZZi * ZY;
w = sqrt(real(conj(w).*w));

indx_planes = [1:n_echo];    % Which planes to consider for fitting
for i = 1 : n_peaks
    % Method 1. Fit a linear function to the logarithms
    linFit = polyfit(t_echo(indx_planes), log(w(i, indx_planes)), 1);
    beta_hat_linFit(i) = -linFit(1);
    % Method 2. Fit an exponential directly
    expFit = fit(t_echo(indx_planes)',w(i, indx_planes)','exp1', 'Exclude', isnan(w(i, indx_planes)));
    beta_hat_expFit(i) = -expFit.b;
end;
beta_hat = beta_hat_expFit(:);

% Comparison figure
R = corrcoef(orig_beta(~isnan(orig_beta)), beta_hat(~isnan(orig_beta)));
figure;
%scatter(orig_beta, beta_hat, 'b.');
plot(orig_beta(1:numel(beta_hat)), beta_hat, 'r^', 'MarkerSize', 3);
hold all;
% scatter(orig_beta, beta_hat_expFit(:), 'b.');
line([0 15], [0 15]);
hold all;
xlabel('FT integration');
ylabel('Model-based');
axis square;
grid on;
% xlim([0.5 3.0]); ylim([0.5 3.0]);    % For T1
% xlim([0 10]); ylim([0 10]);          % For T1rho
% xlim([1 11]); ylim([1 11]);          % For T2
