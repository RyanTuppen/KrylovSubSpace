MinAlpha = 1;%input('What is the Minimum Alpha you want to test? ');
MaxAlpha = 15;%input('What is the Maximum Alpha you want to test? '); %Make it the same as MinAlpha if you want to only use 1 alpha
AlphaSteps = 2;%input('What is the gap between Alphas you want to test (i.e. 1 -> 1.1)? ');
GridSpace = 1;%input('What is your Grid Spacing? ');
MinDim = 10;%input('What is the Minimum Dimensions you want to test? ');
MaxDim = 150;%input('What is the Maximum Dimensions you want to test? '); %Make it the same as MinDim if you want to only have 1 Dimension

I = 1;

NumDim = (MaxDim - MinDim) + 1;%calcs difference
NumAlpha = ((MaxAlpha - MinAlpha)/AlphaSteps) + 1;%calcs alpha num
NumMats = NumDim * NumAlpha;%calcs num mat

Sets = zeros(NumMats,5);%increase the (currently) 5 to 2 + X, where X is however many tests you run

for Dim = MinDim:MaxDim
    for Alpha = MinAlpha:AlphaSteps:MaxAlpha
        A = make_A(Dim,Alpha,GridSpace); %This finds our matrix A
        b = zeros(1,Dim);%sets b

        for i = 1:Dim%loops dimension number of times
            b(i) = sin(i*GridSpace);%sets b to sin function
        end

        Sets(I,1) = Dim; %records the size of the matrix
        Sets(I,2) = Alpha;

        Method1 = tic;%starts time
        GMRES(A,b)
        Sets(I,3) = toc(Method1);%sets time

        Method2 = tic;%starts time
        ConjugateGrad(A,b)
        Sets(I,4) = toc(Method2);%sets time

        Method3 = tic;%starts time
        BiConjGrad(A,b)
        Sets(I,5) = toc(Method3);%sets time
        I = I+1;
    end
end
%Then Save Sets as a .csv - code is:
% writematrix(Sets,'Times.csv')
% can use that in either your preferred spreadsheet program, in R or even Matlab
