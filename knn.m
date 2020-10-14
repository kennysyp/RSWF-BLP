 function [knn_network] = knn(P,k)
    [~,n] = size(P);
    matrix = P;
    A = zeros(n,1);
    knn_network = zeros(n, n);
    P= P-diag(diag(P)); 
   
    [sort_network,idx]=sort(P,2,'descend');
    for i = 1 : n
        A(i)=sum(sort_network(i,1:k));
        knn_network(i,idx(i,1:k))=sort_network(i,1:k);
        %knn_network(i,idx(i,1:k)) = exp(-(norm(matrix(i,:) - matrix(idx(i,1:k),:)).^2) / 4);
        knn_network(i,:) = knn_network(i,:) ./ A(i);
        
    end
 
end