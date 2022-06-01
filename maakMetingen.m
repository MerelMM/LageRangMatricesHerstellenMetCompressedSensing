function [x,meting] = maakMetingen(aantalMetingen,mn, vector)
meting = zeros(1,aantalMetingen);
x=randn(aantalMetingen, mn);
meting = x*vector;
end