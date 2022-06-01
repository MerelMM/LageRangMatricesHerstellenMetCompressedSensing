%MaakMetingenMetSPaarseMatrix
function [x,meting] = MaakMetingenMetPlusMinMatrix(aantalMetingen,mn, vector)
x=maakPlusMinMatrix(aantalMetingen, mn);
meting = x*vector;
end