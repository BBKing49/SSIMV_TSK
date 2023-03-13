function [kernel_width] = kernel_width(C,data,min_kernel,max_kernel)

% U = membership_matrix(C,data);

kernel_width = zeros(size(C,1),size(data,2));
n_samples = size(data,1);
for j=1:size(data,2)
    for k=1:size(C,1)
        kernel_width(k,j) = norm(data(:,j)-ones(n_samples,1)*C(k,j));
    end
end

kernel_width_sum = repmat(sum(kernel_width),size(C,1),1);
kernel_width = kernel_width./kernel_width_sum;

kernel_width = adjustScale(kernel_width,min_kernel,max_kernel);

end

function U = membership_matrix(C,data)

norm_matrix = zeros(size(C,1),size(data,1));
for i=1:size(data,1)
    data_vector = data(i,:);
    for k=1:size(C,1)
        norm_matrix(k,i) = norm(data_vector-C(k,:));
    end
end
norm_sum = repmat(sum(norm_matrix),size(C,1),1);
U = norm_matrix./norm_sum;
end

function outputMatrix = adjustScale(inputMatrix,minVal,maxVal)
minmum=min(min(inputMatrix));
maxmum=max(max(inputMatrix));
outputMatrix=(inputMatrix-minmum)/(maxmum-minmum);
outputMatrix=minVal+outputMatrix*(maxVal-minVal);
end

