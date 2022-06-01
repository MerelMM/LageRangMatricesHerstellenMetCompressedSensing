function[x, iteraties,slechteAB, slechteAB5,slechteAB10] = recover(rijen, kolommen, rang, aantalMetingen)
global iteraties;
%maak originele matrix en metingen
A = randn(rijen,rang);
B = randn(rang, kolommen);
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);

%A<<B
AO = randn(rijen,rang);
BO = randn( rang, kolommen);
A = AO;
B=BO;
A1 = AO/1e1;
B1 = BO*1e1;
z0 = [A1(:);B1(:)];
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-16, 'FunctionTolerance', 1e-16 );
[z] = fminunc(@(z)f(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
slechteAB = norm(matrix-g) /norm(matrix);

A1 = AO/1e2;
B1 = BO*1e2;
z0 = [A1(:);B1(:)];
[z] = fminunc(@(z)f(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
slechteAB5 = norm(matrix-g) /norm(matrix);

A1 = AO/1e3;
B1 = BO*1e3;
z3 = [A1(:);B1(:)];
[z] = fminunc(@(z)f(C,z, kolommen, rijen, rang, meting),z3,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
slechteAB10 = norm(matrix-g) /norm(matrix);
%nu met de oorspronkelijke A en B
iteraties = 0;
z0 = [A(:);B(:)];
[z] = fminunc(@(z)f(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);

function [fun] = f(C,z, kolommen, rijen, rang, meting)
fun= (0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2);
iteraties = iteraties+1;
end
end