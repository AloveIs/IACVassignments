function [ segmentation, centers ] = kmeans_segm(image, K, L, seed)
    width = size(image(:,:,1),2);
    height = size(image(:,:,1),1);
    n_pixels = width*height;
    
    I = double(reshape(image,n_pixels , 3));
    
    
    centers = normrnd(128,64,K,3);
    
    
    for iteration = 1 : L
        segmentation = compute_segmentation(I, centers);

        centers = update_centers(segmentation,I,centers);

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