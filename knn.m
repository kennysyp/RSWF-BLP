 function [knn_network] = knn(P,k)
    [~,n] = size(P);
    matrix = P;
    A = zeros(n,1);
    knn_network = zeros(n, n);
    [sort_network,idx]=sort(P,2,'descend');
    for i = 1 : n
       knn_network(i,idx(i,1:k))=sort_network(i,1:k);
       knn_network(i,:) = knn_network(i,:) ./ A(i);
    end
 
end