function [yT, c_ref, f_ref] = loadData_115x(example)
% A function to load the relaxation experiments data. Takes as its input an
% integer encoding a particular example:
% 1150 - T1
% 1151 - T1rho
% 1155 - T2
%
% Returns:
% yT - the time-domain signal
% c_ref -- spectrometer frequencies in MHz in both dimensions
% f_ref -- frequency offsets in both spectral dimensions

%% Open the Azara file
fid = fopen(fullfile('_data', [num2str(example), 'my_t1t2.spc']));
procdat = fread(fid,inf,'float',0,'l');
fclose(fid);

%% Read the data
dim1 = 1024; %direct dimension points (we get these for "free")
dim2 = 256;  %indirect dimension points (these we would want to undersample)
%dim3 = numel(t_echo);    %time delays for relaxation fit
switch example
    case 1150  % T1
        dim3 = 8;
    case 1151  % T1rho
        dim3 = 10;
    case 1155  % T2
        dim3 = 10;
end;

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
nf1 = 1024;        nf2 = 1024;       % Number of points in each dimension of the spectrum (after zero-filling)
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
%% Cut only an interesting part of the spectrum (i.e. 6.2...10.2 ppm in the 1H direction)
ind_cut = [234:361];      % Indices in the frequency domain to cut (an even number)          ind_cut = [80:207];
yt1 = NaN(nt1, 2*nt2, dim3);     % Data in time domain, complex values along the 1st dimension, R-I-R-I-... along the second
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

% Update the parameters
c_ref = [c1_ref; c2_ref];
f_ref = [f1_ref; f2_ref];

end

