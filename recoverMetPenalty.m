function[x, iteraties,slechteAB, slechteAB5,slechteAB10] = recovermetPenalty(rijen, kolommen, rang, aantalMetingen, lambda)
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
A1 = AO*1e-3;
B1 = BO*1e3;
z0 = [A1(:);B1(:)];
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-16, 'FunctionTolerance', 1e-16 );
[z] = fminunc(@(z)fun(C,z,rijen, rang, kolommen,meting, lambda),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
slechteAB = norm(matrix-g) /norm(matrix);

A2 = AO*1e-4;
B2 = BO*1e4;
z0 = [A2(:);B2(:)];
[z] = fminunc(@(z)fun(C,z,rijen, rang, kolommen,meting, lambda),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
slechteAB5 = norm(matrix-g) /norm(matrix);

A3 = AO*1e-5;
B3 = BO*1e5;
z3 = [A3(:);B3(:)];
[z] = fminunc(@(z)fun(C,z,rijen, rang, kolommen,meting, lambda),z3,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
slechteAB10 = norm(matrix-g) /norm(matrix);

%met oorspronkelijke A en B
iteraties = 0;
z0 = [A(:);B(:)];
[z] = fminunc(@(z)fun(C,z,rijen, rang, kolommen,meting, lambda),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);

function[f] = fun(C,z,rijen, rang, kolommen,meting, lambda)
iteraties = iteraties +1;
f= (0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2) + ...
    lambda * (norm(reshape(z(1:rijen*rang,:), rijen, rang), "fro")^2 - norm(reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen), "fro")^2)^2;
end
end
