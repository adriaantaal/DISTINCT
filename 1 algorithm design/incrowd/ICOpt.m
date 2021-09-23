function [x,f1,f2] = ICOpt(A,y,lambda,x0,max_contrast_ratio,verbose,interiorSolveOption,addSize)
%edited by adriaan taal for ease of use
%  ICOpt: Solution to basis pursuit using the in crowd heuristic
%  x = ICOpt(A,y,[Lambda=1e-5],[max_contrast_ratio=inf],[verbose=0],[interiorSolveOption=1],[addSize=25])
%  Returns the solution "x" to the problem:
%  minimize 1/2 ||A*x - y||^2_2 + Lambda*||x||_1
%  where elements of "x" can be positive or negative.
%
%  Inputs:
%     A: the forward transform, transforming candidate solution
%       "x" into corresponding expected outputs (i.e. y
%       should be close to A*x in the L2 sense for error to be
%       minimized.
%
%     y: the observations you are trying to fit
%
%     Lambda: the regularization parameter.  Defaluts to 1e-5.
%
%     [max_contrast_ratio]: set to 0 any part of x that is less than a
%       factor of max_contrast_ratio smaller than the maximum of the
%       current x solution.  Defaults to infinity.
%
%     [verbose]: controls the verbosity of the code
%
%     [interiorSolveOption]: selects the subproblem solver
%       
%     [addSize]: the number of elements to add every iteration

%0.  Initialize: the "in" crowd is empty.  All the observed signal is
%    unexplained.
s = warning('OFF','optim:quadprog:SwitchToMedScale');
if ~exist('verbose','var')
    verbose = 0;
end
if ~exist('lambda','var')
    lambda = .00001;  %Regularization parameter
end
if ~exist('max_contrast_ratio','var')
    max_contrast_ratio = inf;%  Candidate components of x than this factor less than the largest in magnitude will be set to 0.
end
if ~exist('interiorSolveOption','var')
    interiorSolveOption = 1;
end
status = '';
%Y and Lambda can be converted to double here; with A we can actually live
%with a single-precision matrix and convert to double later.
y = double(y);
lambda = double(lambda);


R = y;

lambdaCut = lambda * (1 + 1e-6);  %  Contributions with a usefullness less than this will not be added to the in crowd; slightly larger than Lambda for numerical reasons.

if ~exist('addSize','var')
    addSize = min([25 min(size(A))]);  %25, unless the problem size is *tiny* - added for compatability.
end


if ~exist('x0','var')
    x = zeros(1,size(A,2));
    in_crowd = [];
    signs = [];
else
    if length(x0) ~= size(A,2)
        if numel(x0)>1
            warning(['Warning from ICOpt:  Starting value of X must be ' num2str(size(A,2)) ' x 1 vector to be compatible with the ' num2str(size(A,1)) ' x ' num2str(size(A,2)) ' A matrix you provided.']);
            disp('Using a zero vector instead.')
        end
        x = zeros(1,size(A,2));
        in_crowd = [];
        signs = [];
    else
        x = x0;
        in_crowd = find(x~=0);
        signs = sign(x(in_crowd));
    end
end

old_in_crowd = 0; %junk start value, != any real value.
old_signs = 0;
run = 1;
for iStep = 1:10000
%     disp(['iStep = ' num2str(iStep)])
    if run == 1

        %%1.  Find L columns of A that could account best
        %%for R.  Add them to the "in" crowd if membership will
        %%not exceed the limit.
        usefulness = R' * A;
        if verbose
            disp(['Max Usefulnes = ' num2str(max(abs(usefulness))) '.']);
        end

        [values,order] = sort(abs(usefulness));
%         addSize = addSize;  %  Number of additions to the "in" crowd; you *could* make this a function of in crowd size if you like.
        %%2. Add only those elements that have a usefulness greater than lambdaCut
        %the subtraction of lambdaCut implies the l1 shrinkage
        comps = values(end:-1:(end-addSize+1)) - lambdaCut;
        shouldAdd = comps > 0;
        toAdd = order(end:-1:(end-addSize+1));
        toAdd = toAdd(shouldAdd);
        in_crowd = [in_crowd toAdd]; %  Adding to the "in" crowd
        if verbose
            disp(['Cardinality of In Crowd = ' num2str(numel(in_crowd)) '.']);
        end

        signs = [signs sign(usefulness(toAdd))];

        %%%%Uncomment the following code if you do not trust your interior
        %%%%solver to get solutions with a gradient of the L2 term less than lambdaCut
        %%%%(see above) on the in crowd set.

        %         if length(unique(in_crowd))<length(in_crowd)
        %             if verbose > 2
        %                 disp(['Notice: duplicate entries in "in" crowd.  Removing...']);
        %             end
        %             in_crowd = unique(in_crowd);
        %         end

        %%3.  Do exact basis pursuit on the "in" crowd.  Some modeled
        %%source magnitudes may become 0.
        Asmaller = A(:,in_crowd);
        Asmaller = Asmaller .* (ones(size(Asmaller,1),1)*signs); %  Flips the sign of the columns of A for which the corresponding entry in "signs" is negative.
        xo = x(in_crowd);
        % Capitol X is the solution on the in crowd, as opposed to x which
        % is the solution on the full space
        % intOpt returns only nonnegative solutions to X.
        if isempty(in_crowd)
            X = zeros(0,1);
        else
            X = intOpt(double(Asmaller),y,lambda,xo,interiorSolveOption);
        end

        %  Set near 0 to 0.
        maxX = max(X);
        if maxX > 0
            X(find(X < maxX/max_contrast_ratio)) = 0;
        else
            X = 0*X; %Set all to 0 if maxX <= 0
        end
        SUPPORT = X~=0;
        %%4.   Kick out all members of the "in" crowd with magnitude 0.
        in_crowd = in_crowd(SUPPORT);
        signs = signs(SUPPORT);

        %%5.   Check for convergence
        if length(in_crowd) == length(old_in_crowd)
            if all((in_crowd.*signs) == (old_in_crowd.*old_signs))
                run=0;
            end
        end


        if iStep>1
            x(old_in_crowd) = 0; %Sets all values of x to 0; since x is sparse this is faster than re-initializing.
        end
        x(in_crowd) = X(SUPPORT);

        if run==0
            break;
        end
        %%6.   Find the R by subtracting the expected
        %%     from the "in" crowd from the observed signal.
        R = y -  Asmaller*X;

        old_in_crowd = in_crowd;
        old_signs = signs;
        
        if nargout > 1
            f1(iStep) = lambda*norm(x,1);
            f2(iStep) = norm(R,'f').^2;
        end

    end
end

%  Add the sign term

x(in_crowd) = x(in_crowd).*signs;
x = x';


if all(x == 0)
    warning(['Warning from "in" crowd optimization:' char(10) ...
        'All guesses are 0.  Probably Lambda is too high, or you''re in the dark.']);
end