function [signal, inst_Frequency, inst_level]=random_band(Fs, N, fc, duration, num_channels, make_plot)
% % random_band: Outputs a test signal with band limited random spectral content.
% %
% % Syntax:  
% %
% % [signal, inst_Frequency, inst_level]=random_band(Fs, N, fc, duration, num_channels, make_plot);
% %
% % **********************************************************************
% %
% % Description
% % 
% % [signal]=random_band;
% % returns a test signal (signal) usings the default input values.  
% % 
% % 
% % The test signal has an instantaneous frequency which randomly
% % varies within the specified frequency range.  The instantaneous level 
% % is constant, the maximum amplitude is 1.   
% % 
% % The test signal is constructed by generating a a set of random number 
% % for the instantaneous frequency as a function of time.  
% % 
% % The instantaneous frequency is integrated to yield the phase angle. 
% % 
% % the trigonometric sin function is used to construct the real part of the 
% % amplitude.  
% % 
% % Filtering and normalizing is used to remove artifacts from the random 
% % number generation.  The resulting test signal is not a warble tone and 
% % it is not pink, white, or brown noise.  
% % 
% % The test signal has characteristics similar to a warble tone.   
% % The level is more constant and the frequuency is band limited
% % like a warble tone.  There is still some random variation in 
% % instantaneous frequency like white noise.  
% % 
% % [signal, inst_Frequency, inst_level]=random_band(Fs, N, fc, duration,
% % num_channels, make_plot);  returns a test signal with a specified
% % samping rate Fs, spectral width of 1/N octave band, center frequency
% % fc, duration in seconds, number of channels num_channels, and will make
% % three plots if make_plot=1;
% %
% % 
% % 
% % See the sections below named Input Variables and Output Variables 
% % for more details
% % 
% % 
% % 
% % **********************************************************************
% %
% % Input Variables
% %
% % Fs=50000;           % (Hz) sampling rate. 
% %                     % default is Fs=50000;
% %
% % N=3;                % Number of bands per octave
% %                     % default is N=3;
% %
% % fc=1000;            % (Hz) center frequency of the test signal
% %                     % forces the user to use the nearest 
% %                     % preferred Nth octave Band frequency.
% %                     % default is fc=1000;  
% %
% % duration=1;         % (seconds) Length of time of the test signal 
% %                     % default is duration=1;
% %
% % num_channels        % Number of channels for the test signal
% %                     % default is num_channels=1;
% %
% % make_plot=1;        % 
% %                     % default is make_plot=1;
% %
% % **********************************************************************
% %
% % Output Variables
% %
% % signal is the time record of the test signal.  The size of the matrix 
% %                      is (num_samples, num_channles)
% %
% % inst_Frequency is the time record of the instantaneous frequency of the 
% %                      test signal.  
% %
% % inst_level is the time record of the instantaneous level of the 
% %                      test signal.   
% % 
% % 
% % **********************************************************************
%
% Example='1';
%
% Fs=44100;         % (Hz) sampling rate
%
% N=3;              % Width of the frequency range is one-third of an 
%                   % octave band 
%
% fc=1000;          % (Hz) center frequency 
%
% duration=5;       % (seconds) of test signal
%
% num_channels=1;   % 1 channel of output for the test signal
%
% make_plot=1;      % plot the instantaneous frequency, 
%                   % instantaneous level, and make a histogram
%                   % of the frequency distribution.  
% 
% [signal, inst_Frequency, inst_level]=random_band(Fs, N, fc, duration, num_channels, make_plot);
%
%
% % **********************************************************************
% %
% % 
% % List of Dependent Subprograms for 
% % random_band
% % 
% % FEX ID# is the File ID on the Matlab Central File Exchange
% % 
% % 
% % Program Name   Author   FEX ID#
% %  1) convert_double		Edward L. Zechmann			
% %  2) envelopes2		Edward L. Zechmann			
% %  3) estimatenoise		John D'Errico		16683	
% %  4) filter_settling_data		Edward L. Zechmann			
% %  5) moving		Aslak Grinsted		8251	
% %  6) nth_freq_band		Edward L. Zechmann			
% %  7) Nth_oct_time_filter2		Edward L. Zechmann			
% %  8) Nth_octdsgn		Edward L. Zechmann			
% %  9) sd_round		Edward L. Zechmann			
% % 10) sub_mean		Edward L. Zechmann			
% % 11) wsmooth		Damien Garcia		NA	
% %
% %
% % **********************************************************************
% %
% % written by Edward L. Zechmann
% %
% % date 4 November 2009
% %
% % modified 13 November    2009
% %
% % modified 23 November    2009    Increase the start and stop band of the 
% %                                 envelope.  Updated Comments.
% %
% % modified 31 December    2009    Automated the number of random jumps 
% %                                 in instantaneous frequency.  
% %                                 The center frequency divided by 10.  
% %                                 
% %
% %
% % **********************************************************************
% %
% %
% % Please Feel Free to Modify This Program
% %
% % See Also: wavwrite, wavread
% %



if (nargin < 1 || isempty(Fs)) || ~isnumeric(Fs)
    Fs=50000;
end

if (nargin < 2 || isempty(N)) || ~isnumeric(N)
    N=3;
end

if (nargin < 3 || isempty(fc)) || ~isnumeric(fc)
    fc=1000;
end

if (nargin < 4 || isempty(duration)) || ~isnumeric(duration)
    duration=5;
end

if (nargin < 5 || isempty(num_channels)) || ~isnumeric(num_channels)
    num_channels=1;
end

if (nargin < 6 || isempty(make_plot)) || ~isnumeric(make_plot)
    make_plot=1;
end




% Only use the preferred third octave band center frequencies
[fc] = nth_freq_band(N, fc, fc);



fc_l=fc*2^(-1/(2*N));

fc_u=fc*2^(1/(2*N));


% Allow the instantaneous frequency to vary up to "jumps_per_second" of different
% values per second. 30 is used as a general number of values per second.
%


% Calculate the instantaneous frequency bins
num_samples_in=max([ceil(fc/10*(duration+1)), 2]);
data_buf=rand(num_samples_in, num_channels);


Fi=fc*2.^((2*data_buf-1)./(2*N));
fu=1:num_channels;
fl=1:num_channels;

for e1=1:num_channels;
    fu(e1)=max([max(Fi(:, e1))*1.1, fc_u]);
    fl(e1)=min([min(Fi(:, e1))/1.1, fc_l]);
end


% Readjust the instantaneous frequency 
for e1=1:num_samples_in;
    for e2=1:num_channels;
        if Fi(e1, e2) >= fc
            Fi(e1, e2)=fc+(fu(e2)-fc)*2/pi*atan((2.2)*(Fi(e1, e2)-fc)/(fu(e2)-fc));
        else
            Fi(e1, e2)=fc+(fc-fl(e2))*2/pi*atan((2.2)*(Fi(e1, e2)-fc)/(fc-fl(e2)));
        end
    end
end


dt_in=duration/(num_samples_in-1);
dt_out=1/Fs;

tin=dt_in*(0:(num_samples_in-1));


num_samples=max([ceil(Fs*duration), ceil(duration*dt_in/dt_out)]);

tout=dt_out*(0:(num_samples-1));

max_t=find(tout < max(tin), 1, 'last' );
tout=tout(1:max_t);

% interpolate the instantaneous frequency bins
Fi2 = interp1(tin, Fi, tout);
clear('Fi', 'tin');

% Calculate the phase angle as a function of time
phi=2*pi*cumsum(Fi2/Fs);

num_samples2=length(Fi2);
clear('Fi2')

% make sure the indices are in the correct order
[num_samples num_channels]=size(phi);

if num_channels > num_samples
    phi=phi';
    [num_samples num_channels]=size(phi);
end


% set the rise time of the envelope
rise_width=0.01;

% set the amount of quiet time for the envelop
quiet_ends=0.01;

start_time=quiet_ends;
stop_time=num_samples2/Fs-start_time;

% make an envelope
[Env]=envelopes2(num_samples2, Fs, start_time, stop_time, rise_width, rise_width);

% Compute the the amplitude from the phase angle
% and apply the envelope
signal=sin(phi).*(ones(num_channels, 1)*(Env'))';

clear('phi');


% Apply the the Nth octave filter if reasonable
num_x_filter=1; 
sensor=1;
settling_time=0.1;
filter_program=0;



if (logical(fc > 4*Fs/num_samples) && logical(Fs/fc > 2)) 
    [fc_out, SP_levels, SP_peak_levels, SP_bands]=Nth_oct_time_filter2(signal, Fs, num_x_filter, N, fc, sensor, settling_time, filter_program);

    % squeeze to the correct number of dimensions
    signal=squeeze(SP_bands);
    clear('SP_bands');
    

    % make sure the indices are in the correct order
    [num_samples num_channels]=size(signal);

    if num_channels > num_samples
        signal=signal';
        [num_samples num_channels]=size(signal);
    end

    % normalize the instanteneous amplitude
    signal=(signal./abs(hilbert(signal)));


    % Apply the envelope for the last time
    signal=signal.*(ones(num_channels, 1)*Env')';

end

% Set the Channel Names for the Plots
channel_names=cell(num_channels, 1);

for e1=1:num_channels;
    channel_names{e1}=['Channel ' num2str(e1)];
end

color_array={'b', 'g', 'r', 'c', 'm', 'k', 'y'};

if nargout > 1
    inst_Frequency=Fs/(2*pi)*gradient(unwrap(angle(hilbert(signal)))')';
    if isequal(make_plot, 1)
        close all;
        figure;
        plot(tout, inst_Frequency);
        maxt=max(tout);
        hold on;
        plot([0 maxt], fc_l*[1 1], 'k');
        plot([0 maxt], fc_u*[1 1], 'k');
        xlabel('Time (seconds)', 'fontsize', 18);
        ylabel('Instantaneous Frequency (Hz)', 'fontsize', 18);
        title([num2str(fc), ' Hz Band Frequency Variation'], 'fontsize', 18);
        ylim(fc*[2^(-1/N) 2^(1/N)]);
        legend(channel_names, 'location', 'NorthEast');
    end
end

if nargout > 2
    inst_level=20*log10(abs(hilbert(signal))./0.00002);
    if isequal(make_plot, 1)
        figure;
        plot(tout, inst_level);
        xlabel('Time (seconds)', 'fontsize', 18);
        ylabel('Instantaneous Level (dB ref. 20 \muPa)', 'fontsize', 18);
        title([num2str(fc), ' Hz Band Level Variation'], 'fontsize', 18);
        ylim(5*round(0.2*20*log10(1/0.00002))+[-10 10]);
        legend(channel_names, 'location', 'NorthEast');
    end
end

if isequal(make_plot, 1)
    
    figure;
    h1=[];
    for e1=1:num_channels;
        buf=inst_Frequency(:, e1);
        buf=buf( (buf > fc_l*2^(-1/(2*N))) & (buf < fc_u*2^(1/(2*N))) );
        hist(buf, 100);
        h2 = findobj(gca,'Type','patch');
        h3=setdiff(h2, h1);
        set(h3,'FaceColor',color_array{mod(e1-1, 7)+1});
        h1=h2;
        hold on;
    end
    
    ylim3=get(gca, 'ylim');
    plot( fc_l*[1 1], ylim3, 'k');
    plot( fc_u*[1 1], ylim3, 'k');
    xlabel('Instantaneous Frequency (Hz)', 'fontsize', 18);
    ylabel('Number of Occurrences', 'fontsize', 18);
    title([num2str(fc), ' Hz Band Frequency Distribution'], 'fontsize', 18);
    legend(channel_names, 'location', 'NorthEast');
end

