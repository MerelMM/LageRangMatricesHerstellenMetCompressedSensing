matrix = zeros(10000,10000);
for i=1:1496612
    if ratings(i,2)<=10000
        matrix(ratings(i,1),ratings(i,2)) = ratings(i,3);
    end
end
spy(matrix)