function prob = mixture_prob(image, K, L, mask)

%Let I be a set of pixels and V be a set of K Gaussian components in 3D (R,G,B).
%  Store all pixels for which mask=1 in a Nx3 matrix
NPX = sum(mask(:));
pixels = zeros(NPX,3);
weights = rand(K,1);
weights = weights/sum(weights);

for i = 1:3
    channel = image(:,:,i);
    pixels(:,i) = channel(mask == 1);
end
%  Randomly initialize the K components using masked pixels

probabilities = zeros(NPX, K);
sigmas = 256 * repmat(eye(3,3),1,1,K);

%mu = 256 * rand(K,3);
mu = normrnd(128,128,K,3);
% 
%     for iteration = 1 : 6
%         segmentation = compute_segmentation(pixels, mu);
% 
%         mu = update_centers(segmentation,pixels,mu);
% 
%     end



%  Iterate L times
for iter = 1 : L
    %     Expectation: Compute probabilities P_ik using masked pixels
    for class = 1 : K
        eig(sigmas(:,:,class))
        
        probabilities(:,class) = weights(class)*mvnpdf(pixels,mu(class,:),sigmas(:,:,class));
    end
    % normalize each row
    normconst = sum(probabilities,2);
    if any(normconst == 0)
        disp("dicaneeee");
        exit();
    end
    probabilities = probabilities ./ normconst;
    %     Maximization: Update weights, means and covariances using masked pixels
    if any(isnan(probabilities))
        disp("dicaneeee");
        exit();
    end
    weights = sum(probabilities,1) / NPX;
    
    for class = 1 : K
        pweights = probabilities(:,class);
        mu(class, :) = sum(pixels.*pweights)/sum(pweights);
        centered = pixels - mu(class, :);
        sigmas(:,:,class) = (centered' * (centered.* pweights))/sum(pweights);
    end
    if any(isnan(sigmas))
        disp("dicaneeee");
        exit();
    end
    fprintf("sigma is %g\n", sigmas(1,1,1));
end

    %  Compute probabilities p(c_i) in Eq.(3) for all pixels I.
    for class = 1 : K
        probabilities(:,class) = weights(class)*mvnpdf(pixels,mu(class,:),sigmas(:,:,class));
    end
    prob = zeros(size(mask));
    prob(mask==1) = sum(probabilities,2);
end

