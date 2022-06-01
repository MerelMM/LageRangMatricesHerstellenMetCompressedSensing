function [x,meting] = MaakMetingenMetSpaarseMatrix2(s, aantalMetingen,mn, vector)
x=maakSpaarseMatrix2(s, aantalMetingen, mn);
meting = x*vector;
end