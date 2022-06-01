function matrix = maakSpaarseMatrix2(spaarsheid, rijen, kolommen)
matrix = zeros(rijen, kolommen);
for i= 1:(rijen*kolommen/spaarsheid)
    matrix(randi(rijen*kolommen)) = 1;
end
end