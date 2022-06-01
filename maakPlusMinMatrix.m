%functie die een mxn matrix genereert van rang k
% alle elementen wel kleiner dan 1; wat is hier een geschikte range voor?
%m,n>=k
function matrix = maakPlusMinMatrix(m,n)
matrix = randi(2,m,n);
for i=1:m
    for j = 1:n
        if(matrix(i,j)==2)
            matrix(i,j)=-1;
        end
    end
end
end


