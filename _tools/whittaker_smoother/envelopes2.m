function [Env]=envelopes2(num_pts, Fs, start_time, stop_time, rise_time_width, fall_time_width)
% % envelopes2: Calculates a cosine envelope for processing signals.
% % 
% % Syntax:
% % 
% % [Env]=envelopes2(time_sec, Fs, tstart, tstop, deltaRise, deltaFall, RiseShape, FallShape);
% %
% % ***********************************************************
% %
% % Description
% % 
% % [Env]=envelopes2;
% % Returns an envelope of length 100000 data points.  The envelope will 
% % have cosine tails at each end which bring the envelop smoothly to 0.  
% % The envelop will begin rising from data point 1 through 50.  
% % Data points 99951 through 100000 bring the envelop smoothly back to 0.   
% % The tails  each have a width of 50 points.  
% % 
% % 
% % [Env]=envelopes2(num_pts);
% % Returns an envelope of length num_pts. The number of data points in 
% % the rise and fall widths is 50.  If num_pts is less than 50 then the 
% % envelope is always 1 and there is no rise to 1 or fall from 1. 
% % 
% % 
% % [Env]=envelopes2(num_pts, Fs, start_time, stop_time);
% % Specifies the sampling rate, start ime (seconds) for the rise from 
% % 0 to 1 and the stop time (seconds) for returning to 0 from 1.  
% % If num_pts is less than 50 then then the envelope is always 1 
% % and there is no rise to 1 or fall from 1. 
% % 
% % 
% % [Env]=envelopes2(num_pts, Fs, start_time, stop_time, rise_time_width, fall_time_width);
% % rise_time_width specifies the time it will take to rise from 0 to 1.   
% % fall_time_width specifies the time it will take to fall from 1 to 0.   
% % 
% % The input variables and output variables sections gives more
% % detail on each variable.  
% % 
% % ***********************************************************
% % 
% % Input Variables
% %
% % num_pts=100000;         % Is the number of data points for the envelope.
% %                         % 
% %                         % default is num_pts=100000;
% % 
% % Fs=50000;               % (Hz) Is the sampling rate.
% %                         % 
% %                         % default is Fs=50000;
% % 
% % start_time=0;           % (seconds) is the time when the signal begins
% %                         % to rise above zero. 
% %                         % 
% %                         % default is start_time=0;        
% % 
% % stop_time=num_pts/Fs;   % (seconds) is the time when the signal 
% %                         % reaches zero.  
% %                         % 
% %                         % default is stop_time=num_pts/Fs;   
% % 
% % rise_time_width=num_pts./Fs*0.0005;
% %                         % seconds) is the amount of time for the signal 
% %                         % to rise from 0 to 1. 
% %                         % 
% %                         % default is rise_time_width=num_pts./Fs*0.0005;         
% % 
% % fall_time_width=num_pts./Fs*0.0005;
% %                         % (seconds) is the amount of time in seconds 
% %                         % for the signal to fall from 1 to 0.  
% %                         % 
% %                         % default is fall_time_width=num_pts./Fs*0.0005;
% %
% % ***********************************************************
% %
% % Output Variables
% % 
% % Env is a column vector of numbers which forms the envelope for signal
% % processing.  
% % 
% % ***********************************************************
% 
% 
% Example='1';
% 
% % The program works without any inputs
% % An envelop with 100000 points is returnd.  
% [Env]=envelopes2;
% figure(1); plot(Env);
% 
% 
% Example='2';
% 
% num_pts=100000;
% [Env]=envelopes2(num_pts);
% figure(2); plot(Env);
% 
% 
% Example='3';
% % an example with rather short rise and fall time widths.  
%
% num_pts=100000;
% Fs=50000;
% start_time=0.0000;
% stop_time=num_pts/Fs;
% rise_time_width=0.0005;
% fall_time_width=0.0005;
% 
% [Env]=envelopes2(num_pts, Fs, start_time, stop_time, rise_time_width, fall_time_width);
% figure(3); plot((1/Fs)*(0:(num_pts-1)), Env);
% 
% 
% 
% Example='4';
% % A delayed start and early stop with slow rise time width and slow fall 
% % time width.
%
% num_pts=100000;
% Fs=50000;
% start_time=0.1;
% stop_time=num_pts/Fs-0.1;
% rise_time_width=0.4;
% fall_time_width=0.4;
% 
% [Env]=envelopes2(num_pts, Fs, start_time, stop_time, rise_time_width, fall_time_width);
% figure(4); plot((1/Fs)*(0:(num_pts-1)), Env);
% 
% % % ***********************************************************
% %
% % envelopes2 is written by Edward L. Zechmann
% %
% % Based on envelopes by William Murphy.
% % 
% %     date 20 September  2008
% % 
% %
% % ***********************************************************
% %
% % Please Feel Free to Modify This Program
% %
% % See Also: envelopes, flattopwin, hann, kaiser, window
% %

if (nargin < 1 || isempty(num_pts)) || ~isnumeric(num_pts)
    num_pts=100000;
end

if (nargin < 2 || isempty(Fs)) || ~isnumeric(Fs)
    Fs=50000;
end

if (nargin < 3 || isempty(start_time)) || ~isnumeric(start_time)
    start_time=0;
end

if (nargin < 4 || isempty(stop_time)) || ~isnumeric(stop_time)
    stop_time=num_pts/Fs;
end

if (nargin < 5 || isempty(rise_time_width)) || ~isnumeric(rise_time_width)
    rise_time_width=num_pts./Fs*0.0005;
end

if (nargin < 6 || isempty(fall_time_width)) || ~isnumeric(fall_time_width)
    fall_time_width=num_pts./Fs*0.0005;
end




dt=(1/Fs);

% Calculate the first index of the rise  
ix_start=ceil(start_time./dt);

if ix_start < 0
    ix_start=0;
end

flag1=0;
if ix_start > num_pts
    flag1=1;
    ix_start=num_pts;
end


% Calculate the last index of the fall 
ix_stop=ceil(stop_time./dt);

if ix_stop > num_pts
    ix_stop=num_pts;
end

flag2=0;
if ix_stop < 0
    flag2=1;
end

% Calculate the width in indices of the rise and fall 
ix_rise_width=ceil(rise_time_width ./dt); 
ix_fall_width=ceil(fall_time_width ./dt); 

if ix_start+ix_rise_width > num_pts
    flag1=1;
end

if ix_start+ix_rise_width > num_pts
    flag1=1;
end

% Calculate the first index of the fall  
ix_stop2=ix_stop-ix_fall_width;

if ix_stop2 < 0
    flag2=1;
end

    
% Initialize the envelope
Env=ones(num_pts, 1);

% Make all points before signal starts equal to zero
Env(1:ix_start)=0;

% Modify Env with the cosine rise 
if ~isequal(flag1, 1)
    uno=ones(ix_rise_width, 1);
    Env(ix_start+(1:ix_rise_width), 1)=Env(ix_start+(1:ix_rise_width), 1).*0.5.*(uno+sin(-pi/2.*uno+pi.*(0:(ix_rise_width-1))'./(ix_rise_width-1)));
end

% Modify Env with the cosine fall
if ~isequal(flag2, 1)
    uno=ones(ix_fall_width, 1);
    Env(ix_stop2+(1:ix_fall_width), 1)=Env(ix_stop2+(1:ix_fall_width), 1).*0.5.*(uno+cos(pi.*(0:(ix_fall_width-1))'./(ix_fall_width-1)));
end

% Make all points after signal stops equal to zero
Env(ix_stop:end)=0;

