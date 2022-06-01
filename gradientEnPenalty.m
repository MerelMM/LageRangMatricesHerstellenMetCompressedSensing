aantalMetingen = 25;
rijen= 6;
kolommen= 5;
rang= 2;
%aantalMetingen= 20;
lambda = 10;
A = randn(rijen,rang);
B = randn(rang, kolommen);
matrix = A*B; %maak matrix van rang rang
vector = matrix(:);
[C,meting] = maakMetingen(aantalMetingen,rijen*kolommen, vector);
A = randn(rijen,rang);
B = randn( rang, kolommen);

z0 = [A(:);B(:)];
options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6 );
objfun = @(z)recoverWithGradientEnPenalty(C,z, kolommen, rijen, rang, meting, lambda)
%[f,g] = objfun(z0);
%g2 = g*0;
%for i = 1:length(g)
%    zh=z0;
%    zh(i) = z0(i)+1e-8;
%    g2(i)=(objfun(zh)-objfun(z0))/1e-8;
%end
%norm(g2-g);
%g2./g;
%g2;
[z] = fminunc(objfun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
x = norm(matrix-g) /norm(matrix)

function[f,gradientt] = recoverWithGradientEnPenalty(C,z, kolommen, rijen, rang, meting, lambda)

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