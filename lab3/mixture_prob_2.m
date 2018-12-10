% https://github.com/txzhao/Image-Analysis-Computer-Vision/tree/master/image-segmentation/script
function [prob, centers] = mixture_prob_2(image, K, L, mask)

% L: number of iterations
% prob: probability of the color of pixels


% Let I be a set of pixels and V be a set of K Gaussian components in 3D (R,G,B)
% Store all pixels for which mask = 1 in a Nx3 matrix
nrows = size(image, 1);
ncols = size(image, 2);
Ivec = im2double(reshape(image, nrows * ncols, 3));
Mvec = reshape(mask, nrows * ncols, 1);
Ima = Ivec(find(Mvec == 1), :);
g = zeros(size(Ima, 1), K);
g1 = zeros(nrows * ncols, K);

% Randomly initialize the K components using masked pixels
[segmentation, centers, ~, ~, ~] = kmeans_mix(Ima, K,L, 1400, true, false);
%[segmentation, centers] = kmeans_segm(Ima, K,L);

cov = cell(K, 1);
cov(:) = {rand * eye(3)};

w = zeros(1, K);
for i = 1 : K
    w(i) = sum(segmentation == i) / size(segmentation, 1);
end

% Iterate L times
for i = 1 : L
%   Expectation: Compute probabilities P_ik using masked pixels
    for k = 1 : K
        mean_k = centers(k, :);
        cov_k = cov{k};
        diff = bsxfun(@minus, Ima, mean_k);
        g(:, k) = 1 / sqrt(det(cov_k) * (2 * pi)^3) * exp(-0.5 *...
            sum((diff * inv(cov_k) .* diff), 2));
    end
    p = bsxfun(@times, g, w);
    norm = sum(p, 2);
    p = bsxfun(@rdivide, p, norm);
    
%   Maximization: Update weights, means and covariances using masked pixels
    w = sum(p, 1) / size(p, 1);
    for k = 1 : K
        centers(k, :) = p(:, k)' * Ima / sum(p(:, k), 1);
        diff = bsxfun(@minus, Ima, centers(k, :));
        cov{k} = (diff' * bsxfun(@times, diff, p(:, k))) / sum(p(:, k), 1);
    end
end

% Compute probabilities p(c_i) in Eq.(3) for all pixels I
for k = 1 : K
    diff = bsxfun(@minus, Ivec, centers(k, :));
    g1(:, k) = 1 / sqrt(det(cov{k}) * (2 * pi)^3) * exp(-0.5 *...
            sum((diff * inv(cov{k}) .* diff), 2));
end
prob = sum(bsxfun(@times, g1, w), 2);
prob = reshape(prob, nrows, ncols, 1);
end


function [segmentation, centers, empty, cen_idx, count] = kmeans_mix(Ivec, K, L, seed, RANDSAMP, DEBUG)
delta = 100 * ones(K, 3);
clus_centers = zeros(K, 3);
clus_centers_new = zeros(K, 3);
idx_old = zeros(1, size(Ivec, 1));
empty = false;
threshold = 0.01;
count = 0;
nthreshold = 2;

% Let X be a set of pixels and V be a set of K cluster centers in 3D (R,G,B).
if RANDSAMP == true
    % Randomly initialize the K cluster centers
%     rng(seed, 'twister');
    clus_centers = rand(K, 3);
    
    % or you can choose cluster centers by yourself
%     clus_centers(1, :) = 1/255 * [255, 155, 0];
%     clus_centers(2, :) = 1/255 * [255, 255, 255];
%     clus_centers(3, :) = rand(1, 3);

    % or you can choose cluster centers from image pixels
%     rng('default');
%     idx = randperm(size(Ivec, 1), K);
%     for i = 1 : K
%         clus_centers(i, :) = Ivec(idx(i), :);
%     end
else
    % implement k-means++ algorithm
    [~, clus_init] = kmeansplus_init(Ivec', K, seed);
    clus_centers = clus_init';
end

% Compute all distances between pixels and cluster centers
D = pdist2(clus_centers, Ivec, 'euclidean');

if DEBUG == false
    % Iterate L times
    for i = 1 : L
        % Assign each pixel to the cluster center for which the distance is minimum
        [~, cen_idx] = min(D);
        
        % Recompute each cluster center by taking the mean of all pixels assigned to it
        for j = 1 : K
            n_idx = find(cen_idx == j);
            if size(n_idx, 2) < nthreshold
                clus_centers(j, :) = rand(1, 3);
%                 empty = 1;
%                 n = randsample(size(Ivec, 1), 1);
%                 clus_centers(j, :) = Ivec(n, :);
            else
                clus_centers(j, :) = double(mean(Ivec(n_idx, :)));
            end
        end
        % Recompute all distances between pixels and cluster centers
        D = pdist2(clus_centers, Ivec, 'euclidean');
    end
else
    % debug
    while (max(delta(:)) > threshold)
        [~, cen_idx] = min(D);
        for j = 1 : K
            n_idx = find(cen_idx == j);
            if isempty(n_idx) == true
                empty = 1;
                n = randsample(size(Ivec, 1), 1);
                clus_centers_new(j, :) = Ivec(n, :);
            else
                clus_centers_new(j, :) = mean(Ivec(n_idx, :));
            end
        end
        delta = abs(clus_centers_new - clus_centers);
%         delta = abs(idx_old - cen_idx);
        count = count + 1;
        D = pdist2(clus_centers_new, Ivec, 'euclidean');
        clus_centers = clus_centers_new;
%         idx_old = cen_idx;
    end
end

[~, cen_idx] = min(D);
segmentation = cen_idx;
centers = clus_centers;
end