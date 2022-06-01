function [x,meting] = MaakMetingenMetSpaarseMatrix(aantalMetingen,mn, vector)
x=maakSpaarseMatrix(aantalMetingen, mn);
meting = x*vector;
end