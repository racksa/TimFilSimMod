clear;

make_symmetric = false;
% If true, the swimmer bodies will be converted to an axisymmetric
% approximation so that the resulting analytic expressions can be used to
% construct axisymmetric bodies in simulations.
% If false, the script will use the original shapes.

images = dir('*_BW.tif');
num_images = length(images);

%% Read in and process the images.
im_mats = cell(1,num_images);

min_im_width = Inf;

for n=1:num_images

    im = imread(images(n).name);

    % Rotate all images to have their largest dimension along the x-axis.
    % MATLAB has some built-in functions for this in the Image Processing
    % toolbox -- no need to re-invent the wheel here.
    if ~islogical(im)
        im = imbinarize(im);
    end
    if size(im,3)>1
        im = squeeze(im(:,:,1));
    end
    im_info = regionprops(im, 'Orientation');
    angle = [im_info.Orientation]; % Sometimes it's a struct rather than a scalar. I think this is if it detects small dots etc. which aren't part of the main shape.
    angle = -angle(1);
    im = imrotate(im, angle);

    % Throw away as much of the black border as we can.
    [rows, cols] = find(im);
    c = [min(cols), max(cols)];
    r = [min(rows), max(rows)];
    im = im(r(1):r(2), c(1):c(2));

    min_im_width = min([min_im_width, r(2)-r(1)]);

    % Make sure the widest part of the shape is in the right-hand side of
    % the image (without loss of generality). This is just to ensure that
    % we don't detect identical larvae swimming in opposite directions as
    % having different shapes.
    [rows_lhs, ~] = find(im(:, 1:round(size(im,2)*0.5) ));
    [rows_rhs, ~] = find(im(:, round(size(im,2)*0.5)+1:end ));
    if (max(rows_rhs) - min(rows_rhs)) < (max(rows_lhs) - min(rows_lhs))
        im = fliplr(im);
    end

    im_mats{n} = im;

end

boundary_mats = cell(1,num_images);
num_angle_bins = 80;
bin_locs = 2*pi*(0:num_angle_bins-1)/num_angle_bins;

for n=1:num_images

    % We've accounted for the orientation of the images, now we need to
    % account for the size. We do so by setting all images to the same width.
    im_mats{n} = imresize(im_mats{n}, [min_im_width, NaN]);

    % There is no more interpolation after this point, so we can discard
    % the interiors safely now too.
    boundary = bwboundaries(im_mats{n}, 'noholes');
    boundary = boundary{1};
    boundary_ind = sub2ind(size(im_mats{n}), boundary(:,1), boundary(:,2));
    temp = zeros(size(im_mats{n}));
    temp(boundary_ind) = 1;
    im_mats{n} = temp;

    % Now we want to construct a radial representation of the boundary.
    % This allows us to guarantee that we unambiguously have boundaries
    % even we seek a lower-dimensional, mode-based representation.
    count = zeros(1, num_angle_bins);
    bins = zeros(size(count));
    for m=1:length(boundary_ind)
        y = round(0.5*size(im_mats{n},1)) - boundary(m,1); % y-index goes in opposite direction to y-coordinate
        x = boundary(m,2) - round(0.5*size(im_mats{n},2));
        theta = mod(atan2(y,x), 2*pi);
        if theta >= (bin_locs(end) + 0.5*bin_locs(2))
            which_bin = 1;
        else
            [~, which_bin] = min(abs(theta - bin_locs));
        end
        count(which_bin) = count(which_bin) + 1;
        bins(which_bin) = bins(which_bin) + norm([x y]);
    end

    boundary_mats{n} = bins./(min_im_width*count);

    if make_symmetric

        % Convert to radial functions which are symmetric about theta = pi.
        % This is to ensure we can make axisymmetric bodies from these curves.
        % Comment this out to decompose the 'original' data instead.
        boundary_mats{n}(2:end) = 0.5*(boundary_mats{n}(2:end) + fliplr(boundary_mats{n}(2:end)));
        boundary_mats{n} = boundary_mats{n}(bin_locs <= pi);

    end

end

%% Find the shape (deformation) modes.

if make_symmetric

    im_data = NaN(num_images, ceil((num_angle_bins-1)/2) + 1);

else

    im_data = NaN(num_images, num_angle_bins);

end

for n=1:num_images

    im_data(n,:) = boundary_mats{n};

end

im_avg = mean(im_data,1);

[U,S,V] = svd(im_data - im_avg, 'econ');
A = U*S; % Coefficients
V = V'; % Modes

if make_symmetric

    % Undo the symmetry-assuming parametrisation.
    V = [V fliplr(V(:, 2:end-(1 - mod(num_angle_bins,2))))];
    im_avg = [im_avg fliplr(im_avg(:, 2:end-(1 - mod(num_angle_bins,2))))];

end

%% Find Fourier coefficients to describe the deformation modes.
max_num_wavenumbers = 1 + ceil(0.5*(num_angle_bins-1));

Vhat = fft(V, [], 2)/num_angle_bins;
Vcos = 2*real(Vhat(:,1:max_num_wavenumbers));
Vsin = -2*imag(Vhat(:,1:max_num_wavenumbers));
Vcos(:,1) = Vcos(:,1)/2;
if max_num_wavenumbers ~= 1 + 0.5*(num_angle_bins-1)
    Vcos(:,max_num_wavenumbers) = Vcos(:,max_num_wavenumbers)/2;
end

im_avg_hat = fft(im_avg)/num_angle_bins;
im_avg_cos = 2*real(im_avg_hat(:,1:max_num_wavenumbers));
im_avg_sin = -2*imag(im_avg_hat(:,1:max_num_wavenumbers));
im_avg_cos(:,1) = im_avg_cos(:,1)/2;
if max_num_wavenumbers ~= 1 + 0.5*(num_angle_bins-1)
    im_avg_cos(:,max_num_wavenumbers) = im_avg_cos(:,max_num_wavenumbers)/2;
end


