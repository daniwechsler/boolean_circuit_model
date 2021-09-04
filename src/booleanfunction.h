/*
 * booleanfunction.h
 *
 * Represents a boolean function with k inputs.
 *
 * f_i := 	A list of booleans of size 2^k.
 *
 * inputs := A list booleans of size k.
 *
 * The value of the function for a given input at position
 * decimal(input) of f_i. The function decimal converts input to
 * its decimal representation as:
 *
 * decimal = input[0]*2^0 + input[1]*2^1 + ... + input[k-1]*2^(k-1)
 *
 *  input
 * 0  1	 2	   index	f_i
 * -------------------------
 * F  F	 F		0		T
 * F  F	 T		1		T
 * F  T	 F		2		T
 * F  T	 T		3		F
 * T  F	 F		4		F
 * T  F	 T		5		F
 * T  T	 F		6		T
 * T  T	 T		7   	T
 *
 */

#ifndef SRC_BOOLEANFUNCTION_H_
#define SRC_BOOLEANFUNCTION_H_

#include <vector>

using namespace std;

class boolean_function {
public:
	boolean_function(vector<bool> f_i);
	virtual ~boolean_function();
	bool compute(bool* inputs);

	vector<bool>  f_i;
	int k;
};

#endif /* SRC_BOOLEANFUNCTION_H_ */
