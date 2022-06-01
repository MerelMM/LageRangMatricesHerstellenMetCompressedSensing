% onderzoek naar het effect beginpunt: A heel klein en B heel groot, vs
% ongeveer gelijk
%product beginpunt is niet uniek
%z = A*B; x geeft resultaat terug met A en B random,  y met A klein en B
%groot
function [x,y] = recover2(rijen, kolommen, rang, aantalMetingen)
A = randn(rijen,rang);
B = randn(rang, kolommen);
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);

%random gekozen A en B
A = randn(rijen,rang);
B = randn( rang, kolommen);
fun= @(z)(0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2);
z0 = [A(:);B(:)];
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-16, 'FunctionTolerance', 1e-16 );
[z] = fminunc(fun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);

%(Ax)(x^-1*B) maak A klein en B groot
A = A*1/10^6;
B = B*10^6;
z0 = [A(:);B(:)];
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-16, 'FunctionTolerance', 1e-16 );
[z] = fminunc(fun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
y = norm(matrix-g) /norm(matrix);
end