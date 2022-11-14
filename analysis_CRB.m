% Analysis of movement: Replicating a plot pipeline
% Cristian Barrera
% 11 / 14 / 22

% plot parameters
TL = [0 5];
nr = 2;
nc = 3;
% Read the needed file
addpath('./bmi-500-2022/');
data = read_trc('./bmi-500-2022/lue-spiral.trc');
nameD = names(data);
% Find the position of the variable needed
index = find(strcmp(nameD, 'L.Finger3.M3'));
% Get the X,Y,Z (assuming first position is X)
marker_name = "L.Finger3.M3";
marker_xyz = d{:,find(names(d) == "L.Finger3.M3") + [0:2]};
% Get the time variable
time = data{:,"Time"};
time_secs = normalize(time,'range',[0, 5]);
t_inds = time>min(TL)&time<max(TL);
t_secs = rem(time(t_inds),1)==0;
%% 1. Plot the X, Y, Z data of the L.Finger3.M3 marker.
figure
subplot(nr,nc,1)
plot(time_secs,marker_xyz(:,1),'b')
hold on
plot(time_secs,marker_xyz(:,2),'r')
hold on
plot(time_secs,marker_xyz(:,3),'y')
hold off
xlabel('seconds')
ylabel('mm')
legend('X','Y','Z')
title('Raw Data')
%% 2. Plot the Y-Z front view.
subplot(nr,nc,2)
plot(marker_xyz(:,2),marker_xyz(:,3),'black')
xlabel('Y')
ylabel('Z')
title('Front View')

%% 3. Filter the original dat two ways:
%  a. Between 2 and 20 Hz for the analysis data.
%  b. Low pass at 2 Hz for the data plotted in panel 3.
% sampling freq fs is the reciprocal of the difference between two points
fs = 1/mean(diff(time));
% cutoff frequencies for the filter
fc_hi = 2;
fc_lo = 20;
order = 6; %6th order filter, high pass
[bh,ah]=butter(order,[fc_hi/(fs/2) fc_lo/(fs/2)]);
marker_filter = filtfilt(bh,ah,marker_xyz);

% fvtool(bh,ah);
[bl,al]=butter(order,fc_hi/(fs/2));
marker_lowfilter = filtfilt(bl,al,marker_xyz);

% plot
subplot(nr,nc,3)
plot(marker_lowfilt(:,2),marker_lowfilt(:,3),'black')
xlabel('Y')
ylabel('Z')
title('Low freq. component')
%% 4. Find and plot the first PC.
[coeff,score,latent] = pca(marker_filter);
first_pca = score(:,1)*coeff(:,1)';
subplot(nr,nc,4)
hold on
plot(marker_filter(:,2),marker_filter(:,3),'k')
plot(first_pca(:,2),first_pca(:,3));
hold off
xlabel('Y')
ylabel('Z')
title('High frequency componet and 1st PC')
%% 5. Use a moving max to find the envelope.
% 6. Use a zero crossing detector to estimate the frequency.
% 7. Find the first projection and plot.

proj = marker_filter*coeff(:,1);

% smooth with a savitsky-golay smoother
proj_smooth = smoothdata(proj,'sgolay');

% count zero crossings
zcd = dsp.ZeroCrossingDetector();
numZeroCross = cast(zcd(proj_smooth(t_inds)),"double");
tremorFrequency = (numZeroCross/2)/max(TL);

% get envelope from 25 sample moving average
env_width = 25;
env = movmax(proj_smooth(t_inds),env_width);

% use the median of the moving maximum as the estimator of the amplitude
amp = median(env);
ttl = round(tremorFrequency,1) + " Hz, " + round(2*amp,1) + " mm amplitude";

% plot
subplot(nr,nc,[5 6])
hold on
plot(time,proj,'k.')
plot(time,proj_smooth,'r')
h1 = refline(0,amp);
h2 = refline(0,-amp);
h1.Color = 0.5*[1 1 1];
h2.Color = 0.5*[1 1 1];
xlim(TL)
title(ttl)
ylabel("mm")
xlabel("seconds")