function[gewoon,gradient,penalty,gradEnPen,itgew,itgrad,itpen, itgradEnPen]=allRecs(rijen,kolommen,rang,aantalMetingen,lambda)
global itgew;
global itgrad;
global itpen;
global itgradEnPen;
itgew=0;
itgrad=0;
itpen=0;
itgradEnPen=0;

%maak originele matrix en metingen
A = randn(rijen,rang);
B = randn(rang, kolommen);
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);

%herstelling met meegegeven penalty
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-7, 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-7 );
A = randn(rijen,rang);
B = randn(rang,kolommen);
z0 = [A(:);B(:)];
z1 = [A(:);B(:)];
z2 = [A(:);B(:)];
[z] = fminunc(@(z)funp(C,z,rijen, rang, kolommen,meting, lambda),z0,options);
g1 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
penalty = norm(matrix-g1) /norm(matrix);

%herstelling met gewone functie
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-7, 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-7 );
[z] = fminunc(@(z)fung(C,z, kolommen, rijen, rang, meting),z1,options);
g2 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
gewoon = norm(matrix-g2) /norm(matrix);

%herstelling met meegegven gradiënt
options = optimoptions(@fminunc,  'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-7,'Algorithm','trust-region','SpecifyObjectiveGradient',true, 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-7 );
[z] = fminunc(@(z)recoverWithGradient(C,z, kolommen, rijen, rang, meting),z2,options);
g3 = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
gradient = norm(matrix-g3) /norm(matrix);

%herstelling met meegegeven gradiënt en penalty
options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-7, 'StepTolerance', 1e-7, 'FunctionTolerance', 1e-7 );
objfun = @(z)recoverWithGradientEnPenalty(C,z, kolommen, rijen, rang, meting, lambda)
[z] = fminunc(objfun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
gradEnPen = norm(matrix-g) /norm(matrix);


%%FUNTIES
%los op door gewoon functie
function [fun] = fung(C,z, kolommen, rijen, rang, meting)
fun= (0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2);
itgew = itgew+1;
end

%los op door meegegeven penalty
function[f] = funp(C,z,rijen, rang, kolommen,meting, lambda)
itpen = itpen +1;
f= (0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2) + ...
    lambda * (norm(reshape(z(1:rijen*rang,:), rijen, rang), "fro")^2 - norm(reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen), "fro")^2)^2;
end

%los op door meegegeven gradiënt
function[f,gradient] = recoverWithGradient(C,z, kolommen, rijen, rang, meting)
itgrad = itgrad+1;
f=(0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2);
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

%los op door meegegeven gradiënt en penalty
function[f,gradientt] = recoverWithGradientEnPenalty(C,z, kolommen, rijen, rang, meting, lambda)
itgradEnPen = itgradEnPen+1;
f = 0.5*norm(C * reshape((reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen)),rijen*kolommen,1)-meting)^2 + ...
    lambda * (norm(reshape(z(1:rijen*rang,:), rijen, rang), "fro")^2 - norm(reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen),"fro")^2)^2;
A = reshape(z(1:rijen*rang,:), rijen, rang);
B  = reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen);
%met A% = eiejT
pos = 1;
gradientt  = zeros(rijen*rang+kolommen*rang,1);
%eerst voor A mxr
for j=1:rang
    for i=1:rijen
        ei = zeros(rijen,1);
        ej = zeros(rang,1);
        ei(i) = 1;
        ej(j)=1;
        E = ei*ej';
          gradientt(pos) = dot(C*reshape((A*B), rijen*kolommen, 1)-meting, C*reshape(E*B, rijen*kolommen, 1)) + ...
            4*lambda*(norm(A,"fro")^2-norm(B,"fro")^2)*trace(A'*E);
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
        gradientt(pos) = dot(C*reshape((A*B), rijen*kolommen, 1)-meting, C*reshape(A*E, rijen*kolommen, 1)) - ...
           4*lambda*(norm(A,"fro")^2-norm(B,"fro")^2)*trace(B'*E);        
       pos = pos+1;
    end
end
end
end

