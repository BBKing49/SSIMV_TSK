function [acc_l,acc_s] = get_Ls(N,view_nums, T)

e=7;
ker = struct('type','gauss','width',1);

acc_x = T;
for i = 1:N
%     ind = find(vec2lab(T) == vec2lab(T(i,:)));
%     acc_s(i,ind) = 1;
    x_vec = acc_x(i,:);
    x_dist = x_vec(ones(N,1),:);
    x_dist = sqrt( sum( (x_dist - acc_x).^2, 2 ) );
    [x_sort,index] = sort(x_dist);
    
    for e_nearest = 1:e
        j = index(e_nearest,:);
        acc_s(i,j) = kernel(ker,acc_x(i,:)',acc_x(j,:)');
    end
    
end
acc_d = diag( sum( acc_s));
acc_l = zeros(size(acc_d));
acc_l = acc_d - acc_s;
end
% end
