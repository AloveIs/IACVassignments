function [ segmentation, centers ] = kmeans_segm(image, K, L, seed)
    width = size(image(:,:,1),2);
    height = size(image(:,:,1),1);
    n_pixels = width*height;
    
    I = double(reshape(image,n_pixels , 3));
    
    centers = mvnrnd(mean(I),50*eye(3),K);
%    disp(centers)
    old_centers = centers;
    
    for iteration = 1 : L
        segmentation = compute_segmentation(I, centers);

        centers = update_centers(segmentation,I,centers);
        centers(isnan(centers)) = 0;
%        disp(norm(centers - old_centers,'fro')/K);

%         if norm(centers - old_centers,'fro')/K < 1
%             fprintf("Stopped at iteration %g\n", iteration);
%             segmentation = reshape(segmentation,height,width);
%             return;
%         end
        old_centers = centers;
    end
    segmentation = reshape(segmentation,height,width);
end

function segmentation = compute_segmentation(I, centers)
    
    [~,segmentation] = pdist2(centers,I,'euclidean','Smallest',1);
    segmentation = segmentation';
end

function centers = update_centers(segmentation,I,centers)
    
    for c = 1: size(centers,1)
       mask = segmentation == c;
       n = sum(mask);
       r = sum(I(mask,1))/n;
       g = sum(I(mask,2))/n;
       b = sum(I(mask,3))/n;
       
       centers(c,:) = [r g b];
    end
end