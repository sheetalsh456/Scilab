// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: harpreet.mertia@gmail.com
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt


function [xopt,resnorm,residual,exitflag,output,lambda] = lsqlin (varargin)
	// Solves a linear quadratic problem.
	//
	//   Calling Sequence
	//   x = lsqlin(C,d,A,b)
	//   x = lsqlin(C,d,A,b,Aeq,beq)
	//   x = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
	//   x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)
	//   x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,param)
	//   [xopt,resnorm,residual,exitflag,output,lambda] = lsqlin( ... )
	//   
	//   Parameters
	//   C : a matrix of doubles, represents the multiplier of the solution x in the expression C*x - d. C is M-by-N, where M is the number of equations, and N is the number of elements of x.
	//   d : a vector of doubles, represents the additive constant term in the expression C*x - d. d is M-by-1, where M is the number of equations.
	//   A : a vector of doubles, represents the linear coefficients in the inequality constraints
	//   b : a vector of doubles, represents the linear coefficients in the inequality constraints
	//   Aeq : a matrix of doubles, represents the linear coefficients in the equality constraints
	//   beq : a vector of doubles, represents the linear coefficients in the equality constraints
	//   LB : a vector of doubles, where n is number of variables, contains lower bounds of the variables.
	//   UB : a vector of doubles, where n is number of variables, contains upper bounds of the variables.
	//   x0 : a vector of doubles, contains initial guess of variables.
	//   param : a list containing the the parameters to be set.
	//   xopt : a vector of doubles, the computed solution of the optimization problem.
	//   fopt : a double, the function value at x.
	//   exitflag : Integer identifying the reason the algorithm terminated.
	//   output : Structure containing information about the optimization.
	//   lambda : Structure containing the Lagrange multipliers at the solution x (separated by constraint type).
	//   
	//   Description
	//   Search the minimum of a constrained linear least square problem specified by :
	//   find the minimum of f(x) such that 
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & 1/2||C*x - d||_2^2  \\
	//    & \text{subject to} & A.x \leq b \\
	//    & & Aeq.x \leq beq \\
	//    & & lb \leq x \leq ub \\
	//    \end{eqnarray}
	//   </latex>
	//   
	//   We are calling IPOpt for solving the linear least square problem, IPOpt is a library written in C++. The code has been written by ​Andreas Wächter and ​Carl Laird.
	//
	// Examples
	// //A simple linear least square example
	// C = [0.9501    0.7620    0.6153    0.4057
	//     0.2311    0.4564    0.7919    0.9354
	//     0.6068    0.0185    0.9218    0.9169
	//     0.4859    0.8214    0.7382    0.4102
	//     0.8912    0.4447    0.1762    0.8936];
	// d = [0.0578
	//     0.3528
	//     0.8131
	//     0.0098
	//     0.1388];
	// A = [0.2027    0.2721    0.7467    0.4659
	//     0.1987    0.1988    0.4450    0.4186
	//     0.6037    0.0152    0.9318    0.8462];
	// b = [0.5251
	//     0.2026
	//     0.6721];
	// [xopt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b)
	//    
	// Examples 
	// C = [0.9501    0.7620    0.6153    0.4057
	//     0.2311    0.4564    0.7919    0.9354
	//     0.6068    0.0185    0.9218    0.9169
	//     0.4859    0.8214    0.7382    0.4102
	//     0.8912    0.4447    0.1762    0.8936];
	// d = [0.0578
	//     0.3528
	//     0.8131
	//     0.0098
	//     0.1388];
	// A =[0.2027    0.2721    0.7467    0.4659
	//     0.1987    0.1988    0.4450    0.4186
	//     0.6037    0.0152    0.9318    0.8462];
	// b =[0.5251
	//     0.2026
	//     0.6721];
	// Aeq = [3 5 7 9];
	// beq = 4;
	// lb = -0.1*ones(4,1);
	// ub = 2*ones(4,1);
	// [xopt,resnorm,residual,exitflag,output,lambda] = lsqlin(C,d,A,b,Aeq,beq,lb,ub)
	//
	// Authors
	// Harpreet Singh


	//To check the number of input and output argument
	[lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 4 | rhs == 5 | rhs == 7 | rhs > 10 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [4 6 8 9 10]"), "lsqlin", rhs);
		error(errmsg)
	end

	C = varargin(1);
	d = varargin(2);
	A = varargin(3);
	b = varargin(4);
	nbVar = size(C,2);

	if ( rhs<5 ) then
		Aeq = []
		beq = []
	else
		Aeq = varargin(5);
		beq = varargin(6);
	end

	if ( rhs<7 ) then
		LB = repmat(-%inf,nbVar,1);
		UB = repmat(%inf,nbVar,1);
	else
		LB = varargin(7);
		UB = varargin(8);
	end


	if ( rhs<9 | size(varargin(9)) ==0 ) then
		x0 = repmat(0,nbVar,1)
	else
		x0 = varargin(9);
	end

	if ( rhs<10 | size(varargin(10)) ==0 ) then
		param = list();
	else
		param =varargin(10);
	end

	if (size(LB,2)==0) then
		LB = repmat(-%inf,nbVar,1);
	end

	if (size(UB,2)==0) then
		UB = repmat(%inf,nbVar,1);
	end

	if (type(param) ~= 15) then
		errmsg = msprintf(gettext("%s: param should be a list "), "lsqlin");
		error(errmsg);
	end


	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "lsqlin");
		error(errmsg);
	end

	options = list(	"MaxIter"     , [3000], ...
					"CpuTime"   , [600] ...
	);

	for i = 1:(size(param))/2

		select param(2*i-1)
			case "MaxIter" then
				options(2*i) = param(2*i);
			case "CpuTime" then
				options(2*i) = param(2*i);
			else
				errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "lsqlin", param(2*i-1));
				error(errmsg)
		end
	end

	nbConInEq = size(A,1);
	nbConEq = size(Aeq,1);

	// Check if the user gives row vector 
	// and Changing it to a column matrix


	if (size(d,2)== [nbVar]) then
		d=d';
	end

	if (size(LB,2)== [nbVar]) then
		LB = LB';
	end

	if (size(UB,2)== [nbVar]) then
		UB = UB';
	end

	if (size(b,2)==nbConInEq) then
		b = b';
	end

	if (size(beq,2)== nbConEq) then
		beq = beq';
	end

	if (size(x0,2)== [nbVar]) then
		x0=x0';
	end

	//Check the size of f which should equal to the number of variable
	if ( size(d,1) ~= size(C,1)) then
		errmsg = msprintf(gettext("%s: The number of rows in C must be equal the number of elements of d"), "lsqlin");
		error(errmsg);
	end

	//Check the size of inequality constraint which should be equal to the number of variables
	if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of columns in A must be the same as the number of elements of d"), "lsqlin");
		error(errmsg);
	end

	//Check the size of equality constraint which should be equal to the number of variables
	if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0 ) then
		errmsg = msprintf(gettext("%s: The number of columns in Aeq must be the same as the number of elements of d"), "lsqlin");
		error(errmsg);
	end

	//Check the size of Lower Bound which should be equal to the number of variables
	if ( size(LB,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "lsqlin");
		error(errmsg);
	end

	//Check the size of Upper Bound which should equal to the number of variables
	if ( size(UB,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "lsqlin");
		error(errmsg);
	end

	//Check the size of constraints of Lower Bound which should equal to the number of constraints
	if ( size(b,1) ~= nbConInEq & size(b,1) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of rows in A must be the same as the number of elementsof b"), "lsqlin");
		error(errmsg);
	end

	//Check the size of constraints of Upper Bound which should equal to the number of constraints
	if ( size(beq,1) ~= nbConEq & size(beq,1) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of rows in Aeq must be the same as the number of elements of beq"), "lsqlin");
		error(errmsg);
	end

	//Check the size of initial of variables which should equal to the number of variables
	if ( size(x0,1) ~= nbVar) then
		warnmsg = msprintf(gettext("%s: Ignoring initial guess of variables as it is not equal to the number of variables"), "lsqlin");
		warning(warnmsg);
	end

	//Check if the user gives a matrix instead of a vector

	if ((size(d,1)~=1)& (size(d,2)~=1)) then
		errmsg = msprintf(gettext("%s: d should be a vector"), "lsqlin");
		error(errmsg); 
	end

	if (size(LB,1)~=1)& (size(LB,2)~=1) then
		errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "lsqlin");
		error(errmsg); 
	end

	if (size(UB,1)~=1)& (size(UB,2)~=1) then
		errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "lsqlin");
		error(errmsg); 
	end

	if (nbConInEq) then
		if ((size(b,1)~=1)& (size(b,2)~=1)) then
			errmsg = msprintf(gettext("%s: Constraint Lower Bound should be a vector"), "lsqlin");
			error(errmsg); 
		end
	end

	if (nbConEq) then
		if (size(beq,1)~=1)& (size(beq,2)~=1) then
			errmsg = msprintf(gettext("%s: Constraint should be a vector"), "lsqlin");
			error(errmsg); 
		end
	end

	for i = 1:nbConInEq
		if (b(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of b can not be negative infinity"), "qpipoptmat");
            error(errmsg); 
        end	
	end
    
	for i = 1:nbConEq
		if (beq(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of beq can not be negative infinity"), "qpipoptmat");
            error(errmsg); 
        end	
	end

	//Converting it into Quadratic Programming Problem

	Q = C'*C;
	p = [-C'*d]';
	op_add = d'*d;
	LB = LB';
	UB = UB';
	x0 = x0';
	conMatrix = [Aeq;A];
	nbCon = size(conMatrix,1);
	conLB = [beq; repmat(-%inf,nbConInEq,1)]';
	conUB = [beq;b]' ; 
	[xopt,fopt,status,iter,Zl,Zu,lmbda] = solveqp(nbVar,nbCon,Q,p,conMatrix,conLB,conUB,LB,UB,x0,options);

	xopt = xopt';
	residual = C*xopt-d;
	resnorm = residual'*residual;
	exitflag = status;
	output = struct("Iterations"      , []);
	output.Iterations = iter;
	lambda = struct("lower"           , [], ..
				    "upper"           , [], ..
				    "constraint"      , []);

	lambda.lower = Zl;
	lambda.upper = Zu;
	lambda.constraint = lmbda;

	select status 
		case 0 then
			printf("\nOptimal Solution Found.\n");
		case 1 then
			printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
		case 2 then
			printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
		case 3 then
			printf("\nStop at Tiny Step\n");
		case 4 then
			printf("\nSolved To Acceptable Level\n");
		case 5 then
			printf("\nConverged to a point of local infeasibility.\n");
		case 6 then
			printf("\nStopping optimization at current point as requested by user.\n");
		case 7 then
			printf("\nFeasible point for square problem found.\n");
		case 8 then 
			printf("\nIterates diverging; problem might be unbounded.\n");
		case 9 then
			printf("\nRestoration Failed!\n");
		case 10 then
			printf("\nError in step computation (regularization becomes too large?)!\n");
		case 12 then
			printf("\nProblem has too few degrees of freedom.\n");
		case 13 then
			printf("\nInvalid option thrown back by IPOpt\n");
		case 14 then
			printf("\nNot enough memory.\n");
		case 15 then
			printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
		else
			printf("\nInvalid status returned. Notify the Toolbox authors\n");
		break;
	end

endfunction
