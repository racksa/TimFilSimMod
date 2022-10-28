% Assumes a .csv file holding nothing but data and all lines contain 4 data: observation number, x-coordinate (microns), y-coordinate (microns) and time (milliseconds). Delete column headings etc. before running this script.
clear;

file_name = 'trace_realunits_3';

%% Work out how many times we have data for, and how many data we have at each time
D = load([file_name '.csv']);
[t_realunits, small_in_big_indices, big_in_small_indices] = unique(D(:,4));
num_times = length(t_realunits);
data_per_time = accumarray(big_in_small_indices, 1);

%% Change the view of the data so that the y-coordinate corresponds to the the normal direction and align the bases
% Plot the original data so we can visualise the starting point
figure;
hold on;
colours = gray(num_times+1);
for n=1:size(D,1)
    plot(D(n,2), D(n,3), '*', 'Color', colours(big_in_small_indices(n), :));
end
title('Original data');
axis equal;
axis off;
hold off;

% Here we make two assumptions:
% (1) the first data point Rebecca has provided at each given time is closest to the filament base;
% (2) the second data point in each filament is displaced from its base
% (approximately) normal to the surface.
bases = D(small_in_big_indices, 2:3);
base_theta = mean(atan2(D(small_in_big_indices+1, 3) - bases(:,2), D(small_in_big_indices+1, 2) - bases(:,1))) - pi/2;
beat_dir = [cos(base_theta); sin(base_theta)];

% Shift all bases to zero and then rotate into the new beat direction
originalD = D; % Save the original form of the data
for n=1:size(D,1)
    D(n, 2:3) = D(n, 2:3) - bases(big_in_small_indices(n), :);
    x = D(n, 2:3)*beat_dir;
    y = norm(D(n, 2:3) - x*beat_dir'); % The cilium isn't allowed to come below the surface.
    D(n, 2:3) = [x, y];
end

% Plot the adjusted data so we can compare it to the starting point
figure;
hold on;
for n=1:size(D,1)
    plot(D(n,2), D(n,3), '*', 'Color', colours(big_in_small_indices(n), :));
end
title('Adjusted data');
axis equal;
axis off;
hold off;

%% Identify the locations at fixed arclength fractions along the filament at each time
% Here, we assume that subsequent data belonging to the same time are
% always at larger arclengths; i.e. Rebecca never went back on her self as
% she recorded data.
num_locs = 20; % Including the base at the origin. For reference, Blake (1972) used num_locs = 10.
sample_points = zeros(num_times, 2*num_locs);
tip_ids = [small_in_big_indices(2:end)-1; size(D,1)]; % We get the tips for free, just like the base.
lengths = zeros(1, num_times);
for n=1:num_times
    sample_points(n, end-1:end) = D(tip_ids(n), 2:3);
    for m=1:data_per_time(n)-1
        lengths(n) = lengths(n) + norm(D(small_in_big_indices(n)+m, 2:3) - D(small_in_big_indices(n)+m-1, 2:3));
    end
    dl = lengths(n)/(num_locs-1);
    dl_mult = 1;
    l_lower = 0;
    l_upper = 0;
    for m=1:data_per_time(n)-1
        l_upper = l_upper + norm(D(small_in_big_indices(n)+m, 2:3) - D(small_in_big_indices(n)+m-1, 2:3));
        if (l_lower < dl_mult*dl)&&(l_upper >= dl_mult*dl)
            fac = (dl_mult*dl - l_lower)/(l_upper - l_lower);
            sample_points(n, 2*(dl_mult+1) - 1 : 2*(dl_mult+1)) = D(small_in_big_indices(n)+m-1, 2:3) + fac*(D(small_in_big_indices(n)+m, 2:3) - D(small_in_big_indices(n)+m-1, 2:3));
            dl_mult = dl_mult + 1;
        end
        l_lower = l_upper;
    end
    sample_points(n,:) = sample_points(n,:)/lengths(n); % Scale to unit length
end

% Plot the sample points
figure;
hold on;
for n=1:num_times
    plot(sample_points(n,1:2:end), sample_points(n,2:2:end), '-o', 'Color', colours(n, :));
end
title('Sample points');
axis equal;
axis off;
hold off;

%% Find the Fourier coefficients for the sample points
temp_fig = figure;
hold on;
plot(sample_points(:,end-1), 'k-');
plot(sample_points(:,end), 'b-');
title('How many data per period?');
estimated_data_per_period = input('Input an (integer) estimate for the number of data per period: '); % Make this estimate based on the tip, which should move the greatest distance. We must use a single value for all data or we will introduce artificial phase differences.
bad_period = true;
while bad_period
    try
        delete(h1);
        delete(h2);
    end
    max_num_wavenumbers = 1 + ceil(0.5*(estimated_data_per_period-1));
    cos_coeffs = zeros(max_num_wavenumbers, 2*num_locs);
    sin_coeffs = zeros(max_num_wavenumbers, 2*num_locs);
    for m=0:1
        dhat = fft(sample_points(1:estimated_data_per_period, end-m))/estimated_data_per_period;
        cos_coeffs(1, end - m) = real(dhat(1));
        for k=2:max_num_wavenumbers
            cos_coeffs(k, end - m) = 2*real(dhat(k));
            sin_coeffs(k, end - m) = -2*imag(dhat(k));
        end
        if max_num_wavenumbers ~= 1 + (estimated_data_per_period-1)/2
            cos_coeffs(max_num_wavenumbers, end - m) = cos_coeffs(max_num_wavenumbers, end - m)/2; % Not part of a conjugate pair. Sine coeff will be zero anyway.
        end
    end
    phi = 2*pi*(0 : num_times-1)/estimated_data_per_period;
    x_reconstruction = cos_coeffs(1, end-1)*ones(1, num_times);
    y_reconstruction = cos_coeffs(1, end)*ones(1, num_times);
    for k=2:max_num_wavenumbers
        x_reconstruction = x_reconstruction + cos_coeffs(k,end-1)*cos((k-1)*phi) + sin_coeffs(k,end-1)*sin((k-1)*phi);
        y_reconstruction = y_reconstruction + cos_coeffs(k,end)*cos((k-1)*phi) + sin_coeffs(k,end)*sin((k-1)*phi);
    end
    figure(temp_fig);
    h1 = plot(x_reconstruction, 'k--');
    h2 = plot(y_reconstruction, 'b--');
    title('Is this a good choice for the number of data per period?');
    input_param = input('Is this a good choice for the number of data per period? Press enter if it is, or enter a new estimate for the number of data per period if not: ');
    if isempty(input_param)
        bad_period = false;
        bad_fit = true;
        phi = 2*pi*(0 : estimated_data_per_period-1)/estimated_data_per_period;
        num_wavenumbers = input(sprintf('Enter a number of modes (an integer between 1 and %i) to retain: ', max_num_wavenumbers)); % I don't think there's a good way of doing this automatically -- even with weighting for low mode numbers, how do we stop the code from finding the zero-error solution of using all modes?
        while bad_fit
            figure(temp_fig);
            hold off;
            plot(sample_points(1:estimated_data_per_period, end-1), 'k-');
            hold on;
            plot(sample_points(1:estimated_data_per_period, end), 'b-');
            x_reconstruction = cos_coeffs(1, end-1)*ones(1, estimated_data_per_period);
            y_reconstruction = cos_coeffs(1, end)*ones(1, estimated_data_per_period);
            for k=2:num_wavenumbers
                x_reconstruction = x_reconstruction + cos_coeffs(k,end-1)*cos((k-1)*phi) + sin_coeffs(k,end-1)*sin((k-1)*phi);
                y_reconstruction = y_reconstruction + cos_coeffs(k,end)*cos((k-1)*phi) + sin_coeffs(k,end)*sin((k-1)*phi);
            end
            h1 = plot(x_reconstruction, 'k--');
            h2 = plot(y_reconstruction, 'b--');
            title('Does this number of modes offer a good fit?');
            input_param = input('Does this number of modes offer a good fit? Press enter if it does, or enter a new choice for the number of modes otherwise: ');
            if isempty(input_param)
                bad_fit = false;
                cos_coeffs = cos_coeffs(1:num_wavenumbers, :); % Don't store coefficients we're never going to use.
                sin_coeffs = sin_coeffs(1:num_wavenumbers, :);
            else
                num_wavenumbers = input_param;
            end
        end
    else
       estimated_data_per_period = input_param;
    end
end
close(temp_fig);

% Having settled on how many wavenumbers to use based on the tip, we can
% find the coefficients for all of the other points. We don't have to do
% anything at the base, where everything is zero by construction.
for n=3:(2*num_locs - 2)
    dhat = fft(sample_points(1:estimated_data_per_period, n))/estimated_data_per_period;
    cos_coeffs(1, n) = real(dhat(1));
    for k=2:num_wavenumbers
        cos_coeffs(k, n) = 2*real(dhat(k));
        sin_coeffs(k, n) = -2*imag(dhat(k));
    end
    if (num_wavenumbers == max_num_wavenumbers)&&(max_num_wavenumbers ~= 1 + (estimated_data_per_period-1)/2)
        cos_coeffs(max_num_wavenumbers, n) = cos_coeffs(max_num_wavenumbers, n)/2;
    end
end

%% Fit polynomials in s to the Fourier coefficients
% The constant term in these polynomials will always be zero to comply with
% the boundary condition at the base exactly. Whilst this means the
% polyfit(...) built-in won't work for our purposes, we can easily
% solve the appropriate least-squares system ourselves using "\".
%
% Furthermore, in both Blake (1972) and Ishikawa et al. (2019) the
% polynomial has degree 3. We assume the same here rather than attempt to
% find the best degree according to some error.
x_cos_coeffs = zeros(num_wavenumbers, 3);
x_sin_coeffs = zeros(num_wavenumbers, 3);
y_cos_coeffs = zeros(num_wavenumbers, 3);
y_sin_coeffs = zeros(num_wavenumbers, 3);
Smat = zeros(num_locs-1, 3);
for n=2:num_locs
    Smat(n-1,1) = (n-1)/(num_locs-1);
end
Smat(:,2) = Smat(:,1).^2;
Smat(:,3) = Smat(:,2).*Smat(:,1);
for n=1:num_wavenumbers
    % cos x-coord
    rhs = cos_coeffs(n, 3:2:end)';
    x_cos_coeffs(n, :) = (Smat\rhs)';
    % cos y-coord
    rhs = cos_coeffs(n, 4:2:end)';
    y_cos_coeffs(n, :) = (Smat\rhs)';
    % sin x-coord
    rhs = sin_coeffs(n, 3:2:end)';
    x_sin_coeffs(n, :) = (Smat\rhs)';
    % sin y-coord
    rhs = sin_coeffs(n, 4:2:end)';
    y_sin_coeffs(n, :) = (Smat\rhs)';
end

% Fitting these functions to the beat pattern means we will, in general,
% have lost the unit length property. Whilst we cannot scale our constant
% coefficients to recover unit length at all times, we can scale them such
% that the average length over a period is 1.
num_space_samples_for_length = 200;
tangent_Smat = ones(3, num_space_samples_for_length);
tangent_Smat(2,:) = 2*(0:num_space_samples_for_length-1)/(num_space_samples_for_length-1);
tangent_Smat(3,:) = 0.75*(tangent_Smat(2,:).^2);

num_time_samples_for_length = 300;
s_array = zeros(1, num_wavenumbers);
c_array = zeros(1, num_wavenumbers);
mean_length = 0;
for n=1:num_time_samples_for_length
    t = (n-1)*2*pi/num_time_samples_for_length;
    for m=0:num_wavenumbers-1
        s_array(m+1) = sin(m*t);
        c_array(m+1) = cos(m*t);
    end
    tx = (c_array*x_cos_coeffs + s_array*x_sin_coeffs)*tangent_Smat;
    ty = (c_array*y_cos_coeffs + s_array*y_sin_coeffs)*tangent_Smat;
    % Trapezium rule
    mean_length = mean_length + 0.5*(norm([tx(1) ty(1)]) + norm([tx(end) ty(end)]));
    for m=2:num_space_samples_for_length-1
        mean_length = mean_length + norm([tx(m) ty(m)]);
    end
end
mean_length = mean_length/(num_time_samples_for_length*(num_space_samples_for_length-1));
x_cos_coeffs = x_cos_coeffs/mean_length;
x_sin_coeffs = x_sin_coeffs/mean_length;
y_cos_coeffs = y_cos_coeffs/mean_length;
y_sin_coeffs = y_sin_coeffs/mean_length;

%% Plot the Fourier-least-squares beat for comparison to the original data.
figure;
hold on;
for n=1:estimated_data_per_period
    t = (n-1)*2*pi/estimated_data_per_period;
    for m=0:num_wavenumbers-1
        s_array(m+1) = sin(m*t);
        c_array(m+1) = cos(m*t);
    end
    x = (c_array*x_cos_coeffs + s_array*x_sin_coeffs)*[zeros(3,1), Smat'];
    y = (c_array*y_cos_coeffs + s_array*y_sin_coeffs)*[zeros(3,1), Smat'];
    plot(x, y, 'o-', 'Color', colours(n,:));
end
title('Reconstructed beat');
axis equal;
axis off;

%% Save the coefficients to file for future use
dlmwrite([file_name '_x_cos_coeffs.dat'], x_cos_coeffs, 'delimiter', ' ', 'precision', '%.6e');
dlmwrite([file_name '_y_cos_coeffs.dat'], y_cos_coeffs, 'delimiter', ' ', 'precision', '%.6e');
dlmwrite([file_name '_x_sin_coeffs.dat'], x_sin_coeffs, 'delimiter', ' ', 'precision', '%.6e');
dlmwrite([file_name '_y_sin_coeffs.dat'], y_sin_coeffs, 'delimiter', ' ', 'precision', '%.6e');