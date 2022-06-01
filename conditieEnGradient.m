function [x, iteraties] = conditieEnGradient(conditieV,rijen, kolommen, rang, aantalMetingen)
global iteraties;
iteraties = 0;
A = randn(rijen,rang);
B = randn(rang, kolommen);
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);
A = randn(rijen,rang);
B = randn( rang, kolommen);
A(:,2)= A(:,2)*conditieV;
z0 = [A(:);B(:)];
options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',5000, 'MaxIterations', 5000,...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8, 'FunctionTolerance', 1e-6);
%'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
%    'OptimalityTolerance', 1e-16, 'StepTolerance', 1e-10, 'FunctionTolerance', 1e-16 
[z] = fminunc(@(z)recoverWGradient(C,z, kolommen, rijen, rang, meting),z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix);
iteraties

function[f,gradient] = recoverWGradient(C,z, kolommen, rijen, rang, meting)
iteraties = iteraties + 1;
f = 0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2;
A = reshape(z(1:rijen*rang,:), rijen, rang);
B  = reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen);
%met A% = eiejT
pos = 1;
gradient  = zeros(rijen*rang+kolommen*rang,1);
%eerst voor A mxr
for j=1:rang
    for i=1:rijen
        ei = zeros(rijen,1);
        ej = zeros(rang,1);
        ei(i) = 1;
        ej(j)=1;
        E = ei*ej';
        gradient(pos) = dot(C*reshape((A*B), rijen*kolommen, 1)-meting, C*reshape(E*B, rijen*kolommen, 1));
        pos = pos+1;
    end
end
%dan voor B rxn
for j=1:kolommen
    for i=1:rang
        ei = zeros(rang,1);
        ej = zeros(kolommen,1);
        ei(i) = 1;
        ej(j)=1;
        E = ei*ej';
        gradient(pos) = dot(C*reshape((A*B), rijen*kolommen, 1)-meting, C*reshape(A*E, rijen*kolommen, 1));
        pos = pos+1;
    end
end
end

end