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