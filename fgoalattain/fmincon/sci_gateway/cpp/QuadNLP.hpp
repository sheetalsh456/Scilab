/*
 * Quadratic Programming Toolbox for Scilab using IPOPT library
 * Authors :
	Sai Kiran
	Keyur Joshi
	Iswarya


 * Optimizing (minimizing) the quadratic objective function having any number of variables and linear constraints.
 *
*/

#ifndef __QuadNLP_HPP__
#define __QuadNLP_HPP__

#include "IpTNLP.hpp"
extern "C"{
#include <sciprint.h>

}
using namespace Ipopt;

class QuadNLP : public TNLP
{
	private:
		Index numVars_;			// Number of variables.
	
		Index numConstr_; 		// Number of constraints.

		const Number *qMatrix_ = NULL;	//qMatrix_ is a pointer to matrix of size numVars X numVars_ 
						// with coefficents of quadratic terms in objective function.

		const Number *lMatrix_ = NULL;//lMatrix_ is a pointer to matrix of size 1*numVars_
						// with coefficents of linear terms in objective function.	
	
		const Number *conMatrix_ = NULL;//conMatrix_ is a pointer to matrix of size numConstr X numVars
						// with coefficients of terms in a each objective in each row.

		const Number *conUB_= NULL;	//conUB_ is a pointer to a matrix of size of 1*numConstr_
						// with upper bounds of all constraints.

		const Number *conLB_ = NULL;	//conLB_ is a pointer to a matrix of size of 1*numConstr_ 
						// with lower bounds of all constraints.

		const Number *varUB_= NULL;	//varUB_ is a pointer to a matrix of size of 1*numVar_ 
						// with upper bounds of all variables.

		const Number *varLB_= NULL;	//varLB_ is a pointer to a matrix of size of 1*numVar_
						// with lower bounds of all variables.

		const Number *varGuess_= NULL;	//varGuess_ is a pointer to a matrix of size of 1*numVar_
						// with initial guess of all variables.
	
		Number *finalX_= NULL;		//finalX_ is a pointer to a matrix of size of 1*numVar_
						// with final value for the primal variables.

		Number *finalZl_= NULL;		//finalZl_ is a pointer to a matrix of size of 1*numVar_
						// with final values for the lower bound multipliers

		Number *finalZu_= NULL;		//finalZu_ is a pointer to a matrix of size of 1*numVar_
						// with final values for the upper bound multipliers

		Number *finalLambda_= NULL;	//finalLambda_ is a pointer to a matrix of size of 1*numConstr_
						// with final values for the upper bound multipliers

		Number finalObjVal_;		//finalObjVal_ is a scalar with the final value of the objective.

		int iter_;			//Number of iteration.

		int status_;			//Solver return status
 
		QuadNLP(const QuadNLP&);
		QuadNLP& operator=(const QuadNLP&);
	public:
		/*
		 * Constructor 
		*/
		QuadNLP(Index nV, Index nC, Number *qM, Number *lM, Number *cM, Number *cUB, Number *cLB, Number *vUB, Number *vLB,Number *vG):
			numVars_(nV),numConstr_(nC),qMatrix_(qM),lMatrix_(lM),conMatrix_(cM),conUB_(cUB),conLB_(cLB),varUB_(vUB),varLB_(vLB),varGuess_(vG),finalX_(0), finalZl_(0), finalZu_(0), finalObjVal_(1e20){	}


		/* Go to :

	http://www.coin-or.org/Ipopt/documentation/node23.html#SECTION00053130000000000000
		For details about these below methods.
		*/
		virtual ~QuadNLP();
		virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
								  Index& nnz_h_lag, IndexStyleEnum& index_style);
		virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
									 Index m, Number* g_l, Number* g_u);
		virtual bool get_starting_point(Index n, bool init_x, Number* x,
										bool init_z, Number* z_L, Number* z_U,
										Index m, bool init_lambda,
										Number* lambda);
		virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);
		virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);
		virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);
		virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
								Index m, Index nele_jac, Index* iRow, Index *jCol,
								Number* values);
		virtual bool eval_h(Index n, const Number* x, bool new_x,
							Number obj_factor, Index m, const Number* lambda,
							bool new_lambda, Index nele_hess, Index* iRow,
							Index* jCol, Number* values);
		virtual void finalize_solution(SolverReturn status,
						   Index n, const Number* x, const Number* z_L, const Number* z_U,
						   Index m, const Number* g, const Number* lambda, Number obj_value,
						   const IpoptData* ip_data,
						   IpoptCalculatedQuantities* ip_cq);
		
		const double * getX();		//Returns a pointer to a matrix of size of 1*numVar 
						// with final value for the primal variables.

		const double * getZu();		//Returns a pointer to a matrix of size of 1*numVars
						// with final values for the upper bound multipliers

		const double * getZl();		//Returns a pointer to a matrix of size of 1*numVars
						// with final values for the upper bound multipliers

		const double * getLambda();	//Returns a pointer to a matrix of size of 1*numConstr
						// with final values for the constraint multipliers


		double getObjVal();		//Returns the output of the final value of the objective.

		double iterCount();		//Returns the iteration count

		int returnStatus();		//Returns the status count

		
};

#endif __QuadNLP_HPP__
