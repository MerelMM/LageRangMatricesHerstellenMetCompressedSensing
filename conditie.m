function[x,y ] = recover(rijen, kolommen, rang, aantalMetingen)
global iteraties;
%maak originele matrix en metingen
A = randn(rijen,rang);
A = orth(A);
B = randn(rang, kolommen);
B = orth(B')'
A = A*diag([100,1,1-1e-8]); %SV => 100,1,0 conditiegetal 1e-8
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-16, 'FunctionTolerance', 1e-16 );
rang = rang-1;
A = randn(rijen, rang);
A1 = A
B = randn(rang, kolommen);
z0 = [A1(:);B(:)];
[z] = fminunc(@(z)f(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);

%nu met de oorspronkelijke A en B
iteraties = 0;
z0 = [A(:);B(:)];
[z] = fminunc(@(z)f(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
y = norm(matrix-g) /norm(matrix);

function [fun] = f(C,z, kolommen, rijen, rang, meting)
fun= (0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2);
iteraties = iteraties+1;
end
end