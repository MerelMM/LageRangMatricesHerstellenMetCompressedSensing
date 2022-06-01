%vergelijk verschillende metingen

%afmetingen matrices
drempel = 10e-5;
rijen = 7;
kolommen = 6;
rang = 2;
start = (rijen+kolommen-rang)*rang;
aantalMetingen = 1.5*start;
aantalKeer = 5;
y = zeros(aantalKeer*(aantalMetingen-start+1),1);
k=0;
m=0;
lambda=10;

xRandom = zeros(aantalKeer*(aantalMetingen-start+1),1);
xSparse = zeros(aantalKeer*(aantalMetingen-start+1),1);
xPlusMin = zeros(aantalKeer*(aantalMetingen-start+1),1);

for i= start:aantalMetingen
    for j=1:aantalKeer
        k=k+1;   
        [xRandom(k), xSparse(k), xPlusMin(k)] = losOpMetGradEnPen(rijen,kolommen,rang,i,lambda);
    end
end

g=0;
k=0;
succesRecoveryRandom = zeros((aantalMetingen-start+1),1);
succesRecoverySparse = zeros((aantalMetingen-start+1),1);
succesRecoveryPlusMin = zeros((aantalMetingen-start+1),1);
for j=start:aantalMetingen
    g=g+1;
    random=0;
    sparse=0;
    plusmin=0;
    for i = 1:aantalKeer
        k=k+1;
        if (xRandom(k)<=drempel)
            random = random + 1;
        end
        if (xSparse(k)<=drempel)
            sparse = sparse + 1;
        end
        if (xPlusMin(k)<=drempel)
            plusmin = plusmin + 1;
        end
    end
    succesRecoveryRandom(g) = random/aantalKeer;
    succesRecoverySparse(g) = sparse/aantalKeer;
    succesRecoveryPlusMin(g) = plusmin/aantalKeer;

    % resultaten plotten
    x=start/start:1/start:aantalMetingen/start;
    fig1 = figure(1)
    plot(x,succesRecoveryRandom,'--ks')
    title('metingsmatrix random gegenereerd',fontsize=16)
    xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
    ylabel('succesvol hersteld [%]',fontsize=16)
    grid on

    fig2 = figure(2)
    plot(x,succesRecoverySparse,'--ks')
    title('ijle metingsmatrix',fontsize=16)
    xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
    ylabel('succesvol hersteld [%]',fontsize=16)
    grid on
    
    fig3 = figure(3)
    plot(x,succesRecoveryPlusMin,'--ks')
    title('metingsmatrix met random +1/-1 elementen',fontsize=16)
    xlabel('oversamplingsfactor \phi [aantal metingen / (m+n-rang)*rang]',fontsize=16)
    ylabel('succesvol hersteld [%]',fontsize=16)
    grid on
end

%%FUNTIES
function [xRandom, xSparse, xPlusMin] = losOpMetGradEnPen(rijen,kolommen,rang,aantalMetingen,lambda)
%maak matrix
A= randn(rijen, rang);
B= randn(rang, kolommen);
matrix = A*B;
vector = matrix(:);

%neem verschillende soorten metingen
[CRandom,metingRandom] = maakMetingen(aantalMetingen,rijen*kolommen, vector);
[CSparse, metingSparse] = MaakMetingenMetSpaarseMatrix(aantalMetingen,rijen*kolommen, vector);
[CPlusMin, metingPlusMin] = MaakMetingenMetPlusMinMatrix(aantalMetingen,rijen*kolommen, vector);

%beginpunt
A = randn(rijen,rang);
B = randn(rang,kolommen);
z0 = [A(:);B(:)];

options = optimoptions(@fminunc,'Algorithm','trust-region','SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations',999999, 'MaxIterations', 1e6,...
    'OptimalityTolerance', 1e-17, 'StepTolerance', 1e-17, 'FunctionTolerance', 1e-17 );

%oplossen met random matrix metingen
objfun = @(z)recoverWithGradientEnPenalty(CRandom,z, kolommen, rijen, rang, metingRandom, lambda);
[z] = fminunc(objfun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
xRandom = norm(matrix-g) /norm(matrix);

%oplossen met sparse matrix metingen
objfun = @(z)recoverWithGradientEnPenalty(CSparse,z, kolommen, rijen, rang, metingSparse, lambda);
[z] = fminunc(objfun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
xSparse = norm(matrix-g) /norm(matrix);

%oplossen met plusminmatrix metingen
objfun = @(z)recoverWithGradientEnPenalty(CPlusMin,z, kolommen, rijen, rang, metingPlusMin, lambda);
[z] = fminunc(objfun,z0,options);
g = (reshape(z(1:rijen*rang,:), rijen, rang) ...
    *reshape(z(rijen*rang+1:(rijen*rang+kolommen*rang),:),rang,kolommen));
xPlusMin = norm(matrix-g) /norm(matrix);

end

%los op door meegegeven gradiënt en penalty
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