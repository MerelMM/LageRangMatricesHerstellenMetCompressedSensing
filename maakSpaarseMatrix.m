function matrix = maakSpaarseMatrix(rijen, kolommen)
matrix = zeros(rijen, kolommen);
for i= 1:(rijen*kolommen/2)
    matrix(randi(rijen*kolommen)) = 1;
end
end
