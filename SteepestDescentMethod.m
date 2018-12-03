
clc
clear all

% Getting inputs from User, namely:
%1. The training data file, 2. Number of Inputs, 3.Number of outputs.
%If the input file is a .tra format file then use the method below to read the
%file:

File = input('Enter the training file name: ', 's');
N = input('Enter the number of inputs (N): ');
M = input('Enter the number of outputs (M): ');

% The following code reads the file and stores all the patterns in an Nv by (N+M) matrix
fid = fopen(File, 'r');
File_Values = fscanf(fid, '%f');
fclose(fid);
Nv = numel(File_Values)/(N+M);
fprintf('Number of patterns (Nv) = %d\n', Nv);
File_Values = reshape(File_Values, [(N+M) Nv])';

% Store the inputs in variable x where x is the Input of dimension [Nv x N]
% and the outputs in variable t where t is the Desired (Target) output of dimension [Nv x M]
x = File_Values(:, 1:N);
t = File_Values(:, N+1:N+M);
clear File_values;



%Alternatively, If the input file is a .mat format file where inputs have been already
%stored in variable x and the target outputs in variable t, then you can use the following steps:
%File = input('Enter the training file name: ', 's');
%load(File);
%N = size(x,2);
%M = size(t,2);
%Nv= size(x,1);



%Now we have values of N, M, Nv  and variables x and t with us.

%For better results, we insert an extra input=1 as the last column of the input vector, i.e N+1th column.

xAug = [x ones(Nv,1)];

%R is the  Autocorrelation matrix of size(N+1) by (N+1)

R = (xAug' * xAug) / Nv;

%C is the Cross-correlation matrix of size(N+1) by (M)

C = (xAug' * t) / Nv;

% Calculates output energies
Et = sum(t .* t) / Nv;


MSE=0;
for k=1:M
    W = zeros(N+1,1); %W= Weights matrix for kth output
    g = zeros(N+1,1); %g= Gradient matrix for kth output
    
    % Loop over iterations for better training
    for iter = 1:100
        
        % Calculation for gradient g
        g = -2 * (C(:,k) - R * W);
        
        % Calculate energy in the gradient
        Eg = sum(g .* g);
        
        % Calculate B2, where B2 is the Optimal Learning Factor
        Num = sum(g .* g)/2;
        Den = g' * R *  g;
   
        B2=Num/Den; 
        
        % Update weights
        W = W - B2*g;
            
    end
    
    %Calculating error at each node.
     E(k) = Et(k) - 2 * W' * C(:,k) + W' * R * W; 
     
     WW(:,k)=W; %WW is the weight matrix for all outputs combined.
end

 fprintf("\nThe weights for the given network are:")
 W=WW
 
 %Adding up the errors to obtain the network's total Mean Square Error.
for m=1:M
    fprintf("Error at node %d = %f\n",m,E(m))
  MSE=MSE+E(m);
end
fprintf("\nMean Square Error of the network =%f\n", MSE)









