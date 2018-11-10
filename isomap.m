numberofdatax=50;
numberofdatay=50;
numberofneighbors=30;
reduceddim=2;
indata = zeros(numberofdatax*numberofdatay, 3);
V = size(indata, 1);
% generate Data Matrix indata, you can replace this part with your own data. 
% noise in the data
noiseScaling = 0.3; 

for i=1:numberofdatax
    for j=1:numberofdatay
        indata(i*(numberofdatay-1)+j,1) = (i/numberofdatax * 2.0 * 3.14 + 3.14) * cos(i/numberofdatax * 2.0 * 3.14) + (rand)*noiseScaling;
        indata(i*(numberofdatay-1)+j,2) = (i/numberofdatax * 2.0 * 3.14 + 3.14) * sin(i/numberofdatax * 2.0 * 3.14) + (rand)*noiseScaling;
        indata(i*(numberofdatay-1)+j,3) = 10.0*j/numberofdatay + (rand)*noiseScaling;
    end
end

% Visualization of the input dataset. Work only for 2D& 3D 
if size(indata, 2) == 3
   scatter3(indata(:,1), indata(:,2), indata(:,3))
elseif size(indata, 2) == 2
   scatter(indata(:,1), indata(:,2))
else
   warning('dataset cannot be visualized')
end 


D = zeros(V, V); 

for i=1:V
    for j = i+1:V
        distance  = norm(indata(i, :)-indata(j, :));
        D(i, j) = distance;
        D(j, i) = distance;
    end
end


%Assign distance 1000 to the non-neighbors(adjust this value if it does not fit in you dataset)

Dsort = sort(D);
for k=1:V
    threshold = Dsort(k, numberofneighbors);
    for l=k+1:V
        if D(k,l)>threshold
           D(k,l)=1000;
           D(l,k)=1000;
        end
    end
end    
 
% Determine the shortest path between all pairs using Floyd-Warshall algorithm 
dataMatrix = FloydWarshall(D);

% Centering J matrix to use in centering
mOne = ones(V, V);
mOne= mOne/V;
J = eye(V, V)-mOne;


% Double centering B matrix
B = -0.5*J*dataMatrix*J;

% Eigendecomposition of B
e = eig(B);
[W, DD] = eig(B);
W = real(W);

% reduced eigenvector matrix
Edim = W(1:reduceddim,:);
deltaSqrt = zeros(reduceddim, reduceddim);

for l=1:reduceddim
    deltaSqrt(l,l)=sqrt(e(l));
end    
    
% Dim-reduced(Isomap based) data
reduceddata = Edim'*deltaSqrt;
    
%scatter(reduceddata(:,1), reduceddata(:,2))