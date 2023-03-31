function [] = BiConjGrad(A,b)%Bi-conjugate method function

	res_tol  = 1e-9;%sets resolution tolerance
	max = 1000; %sets maximum number of iterations
	N = size(A, 1);%Finds size of A
	x = zeros(N, 1);%sets vector x

	r = b - A * x;%calculates r
	p = r;%sets values
	rs = r;
	ps = p;
	r2 = r' * rs;%sets r2 so can be interchanged later
	
	residual = norm(r, 2);%calculates residual norm
	res_stop = residual * res_tol;%stop condition
	count = 1;%initialises count
	res_norm(count) = residual;%sets reolution norm
	
	converged = 0;%initialises converged variable
	while ((count < max) && (residual > res_stop))%loops until stop conditions met
		alpha = r2 / ((A*p)' * ps);%Calculates alpha
		x = x + alpha * p;%adds onto x
		r = r - alpha * (A*p);%Finds r
		rs = rs - alpha * (A' * ps);%Finds rs

		%sets new and old r2
		r2_old = r2;
		r2 = r' * rs;
		beta = r2 / r2_old;
		
        %finds p and ps
		p  = r  + beta * p;
		ps = rs + beta * ps;
		
		count = count + 1;%increments count
		residual = norm(r, 2);%calculates new residual
		res_norm(count) = residual;
	end

    if residual <= res_stop %sets converged to 1 if residual less than stop
        converged = 1; 
   end
end


