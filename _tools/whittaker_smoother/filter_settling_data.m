function [y2, num_settle_pts, settling_time]=filter_settling_data(Fs, y, settling_time)
% % filter_settling_data: Creates data to preappend to a time record for settling a filter
% %
% % Syntax:
% %
% % [y2, num_settle_pts]=filter_settling_data(Fs, y, settling_time);
% %
% % **********************************************************************
% %
% % Description
% %
% % [y2, num_settle_pts]=filter_settling_data(Fs, y, settling_time);
% % Returns y2 a set of data for settling a filter having num_settle_pts of
% % data points.  Fs (Hz) is the sampling rate, y is the time record, and
% % settling_time (s)is the time require for the impulse response of the
% % filter to decay to the desired amount.
% %
% % The filter settling data is the reflection of the data about the first
% % data point.  The settling data is phase matched to the original data
% % and smoothed to remove transient effects.  A sinusoidal input will
% % produce a pure sinusoid with continuous phase when the the filter
% % settling data is preappended to the original signal.
% %
% %
% %
% % ***********************************************************
% %
% % Input Variables
% %
% % Fs is the sampling rate in Hz.  default is 50000 Hz.
% %
% % y is the multichannel input time record in (Pa).  Processsing assumes
% %      that y has more channels than time record samples.
% %      y=randn(10000, 10) is the default.
% %
% % settling_time is the time it takes the filter to settle (seconds).
% %      default is settling_time=0.1;
% %
% % ***********************************************************
% %
% % Output Variables
% %
% % y2 is the filter settling data
% %
% % num_settle_pts is the number of data points in the settling data.
% %
% % **********************************************************************
%
%
% % These examples show the need for using the settling data especially
% % when using the filter command.
%
%
% Example='1';
% % A sinusoid is processed and the change from the filter settling data to
% % the actual signal is almost perfectly continuous.
%
% f_sig=20;
% Fs=50000;
% t_SP=(0:(1/Fs):(100/f_sig));
% y=sin(2*pi*t_SP*f_sig);
% settling_time=0.1;
% [y2, num_pts_se]=filter_settling_data(Fs, y, settling_time);
% buf=[y2 y];
% buf2=hilbert(buf);
% t_SP=1/Fs*(1:length(buf));
% subplot(3,1,1); plot(t_SP, buf); ylabel('signal (Pa)');
% hold on;
% subplot(3,1,2); plot(t_SP, abs(buf2)); ylabel('Amplitude (Pa)');
% subplot(3,1,3); plot(t_SP, 180/pi*angle(buf2)); ylabel('Phase deg.');
%
%
% Example='2';
% The lack of settling the filter causes the beginning of the filtered time
% record to oscillate with a very high amplitude.
%
% load Example_Data
% t=1/Fs*(1:length(SP));
% settling_time=0.1;
%
% [y2, num_pts_se]=filter_settling_data(Fs, SP, settling_time);
%
% % Now apply a 1/3 octave band filter
% Fc=100;
% N=3;
% n=3;
%
% [Bc, Ac]=Nth_octdsgn(Fs, Fc, N, n);
% SP2 = filter(Bc, Ac, [y2 SP]);
% SP2=SP2((num_pts_se+1):end);
% SP22 = filtfilt(Bc, Ac, [y2 SP]);
% SP22=SP22((num_pts_se+1):end);
%
% SP1 = filter(Bc, Ac, [SP]);
% SP11 = filtfilt(Bc, Ac, [SP]);
% plot(t, SP2, 'k');
% hold on;
% plot(t, 0.05+SP22, 'b');
% plot(t, 0.1+SP1, 'r');
% plot(t, 0.15+SP11, 'c');
% legend('"filter" With Filter settling', '"filtfilt" With Filter settling', '"filter" No Filter Settling',  '"filtfilt" No Filter Settling');
%
%
% Example='3';
%
% % Only the first data point using the resample program is problematic,
% % however this discontinuity can cause problems with additional signal
% % processing.
%
% load Example_Data
% [y2, num_pts_se]=filter_settling_data(Fs, SP, settling_time);
% SP2=resample([y2 SP], 1, 8);
% SP2=SP2((floor(num_pts_se/8)+1):end);
% SP1=resample(SP, 1, 8);
% plot(SP2, 'k');
% hold on;
% plot(SP1, 'r');
% legend('With Filter settling', 'No Filter Settling');
%
%
%
%
% % ***********************************************************
% %
% % Subprograms
% % 
% % List of Dependent Subprograms for 
% % filter_settling_data
% % 
% % FEX ID# is the File ID on the Matlab Central File Exchange
% % 
% % 
% % Program Name   Author   FEX ID#
% % 1) convert_double		Edward L. Zechmann			
% % 2) estimatenoise		John D'Errico		16683	
% % 3) wsmooth		Damien Garcia		NA	
% %
% % **********************************************************************
% %
% % Program was written by Edward L. Zechmann
% %
% %      date 6 December    2008
% %
% % modified  8 December    2008    Added Comments and an example.
% %
% % modified 17 December    2008    Added descriptions of the subprograms,
% %                                 input variables, and output variables.
% % 
% % modified  6 October     2009    Updated comments
% % 
% %
% % **********************************************************************
% %
% % Please feel free to modify this code.
% %
% % See Also: filter, filtfilt, resample, ACweight_time_filter,
% %           hand_arm_time_fil, whole_body_time_filter


if (nargin < 1 || isempty(Fs)) || ~isnumeric(Fs)
    Fs=50000;
end

if (nargin < 2 || isempty(y)) || ~isnumeric(y)
    y=rand(1, 50000);
end

% Make the data have the correct data type and size
[y]=convert_double(y);

[m1,n1]=size(y);
flag1=0;

if m1 > n1
    y=y';
    [m1, n1]=size(y);
    flag1=1;
end


if (nargin < 3 || isempty(settling_time)) || ~isnumeric(settling_time)
    settling_time=0.1;
end



% calculate the number of data points needed to settle the filter
num_settle_pts=ceil(Fs*settling_time);


if num_settle_pts < 1 && logical(n1 > 1)

    y2=[];
    num_settle_pts=0;

else

    % initialize the output matrix
    %y2=zeros(m1, num_settle_pts);

    % Calculate the time array
    t=1/Fs*(0:(num_settle_pts-1));

    % Match the amplitude of the signal to the amplitude of a sinusoid.
    A=sqrt(2)*sqrt(sum(y.^2, 2))./sqrt(n1);
    A(A==0)=1;
    y0=y(:,1);

    % Calculate the initial slope
    grad_pts=min([ n1, num_settle_pts]);
    bb=gradient(y(:, 1:grad_pts))*Fs;
    sl=bb(:, 1);
    max_sl=max(abs(bb),[], 2);
    
    % Calculate the signal frequency to be one-hundreth of the sampling rate
    Fsig=max([1/settling_time, sl'./(2*pi)]);
    Fsig_max=max(max_sl./2*pi);

    if Fsig_max > Fs/20

        % initialize the output matrix
        y2=zeros(m1, num_settle_pts);

        % Calculate the number of iterations needed to splice the settling time
        % data together.
        num_iter=ceil(num_settle_pts./n1);

        if num_iter < 1
            num_iter=1;
        end

        max_pts=min([num_settle_pts, n1]);

        w=hann(num_settle_pts);

        for e1=1:m1;

            y1=y(e1, 1:max_pts);

            for e2=1:num_iter;

                num_last_pts=num_settle_pts-(num_iter-1)*max_pts;

                if  isequal(e2, num_iter);
                    nn=num_last_pts;
                else
                    nn=max_pts;
                end


                if isequal(mod(e2, 2), 1)
                    y2(e1, (e2-1)*max_pts+(1:nn))=y1(1:nn);
                else
                    y2(e1, (e2-1)*max_pts+(1:nn))=fliplr(y1((max_pts-nn+1):max_pts));
                end

            end

            y2(e1, :) = w'.*wsmooth(y2(e1, :), 3);

        end

        y2=fliplr(y2);

    else
        
        % Calculate the offet normalized slope
        B=2*pi*Fsig;
        norm_sl=sl./A./B;

        % If the normalized slope is too steep the phase cannot be matched exactly
        norm_sl(norm_sl > 1)=0.99;
        norm_sl(norm_sl < -1)=-0.99;

        % Offset phase angle
        ph0=acos(norm_sl);

        % Calculate the intial phase angle to match the slope at the
        % last point to the actual slope with teh normalized slope being
        % in the range of -1 < norm_sl < 1
        ph=ph0-B.*t((num_settle_pts-1));

        y2=A.*sin(B*t+ph);

        % Shift sinusoid vertically to match y(o) constraint exactly
        y2=y2-y2(:, end)+y0;
    end

end



% Make sure the array has the same orientation as before
% transpose if necessary.
if isequal(flag1, 1)
    y2=y2';
end


