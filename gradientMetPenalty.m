%maakSlechteConditie op derde kolom
function[x, iteraties] = gradientMetConditie(conditieV, rijen, kolommen, rang, aantalMetingen)
global iteraties;
iteraties = 0;
A = randn(rijen,rang);
B = randn(rang, kolommen);
A(:,2)= A(:,2)*conditieV;
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);
A = randn(rijen,rang);

z0 = [A(:);B(:)];
options = optimoptions(@fminunc,'MaxFunctionEvaluations',5000, 'MaxIterations', 5000,...
    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-10, 'FunctionTolerance', 1e-10 );
[z] = fminunc(@(z)recoverHier(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);

function[f] = recoverHier(C,z, kolommen, rijen, rang, meting)
iteraties= iteraties + 1
f=(0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2) + ...
    10 * (norm(reshape(z(1:rijen*rang,:), rijen, rang), "fro")^2 - norm(reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen), "fro")^2)^2;
A = reshape(z(1:rijen*rang,:), rijen, rang);
B  = reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen);

end
end