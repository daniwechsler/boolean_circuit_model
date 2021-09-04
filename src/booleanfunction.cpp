/*
 * booleanfunction.cpp
 */

#include "booleanfunction.h"
#include "circuitmodel.h"
#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

boolean_function::boolean_function(vector<bool> f_i) {

	this->f_i = f_i;
	int size = f_i.size();
	int k = 0;
	while(true) {

		if (size == 1) {
			break;
		}
		if (size % 2 == 0) {
			k += 1;
			size = size / 2;
		} else {
			assert(false);
		}
	}
	this->k = k;
}


boolean_function::~boolean_function() {
}


bool boolean_function::compute(bool* inputs) {

	return this->f_i[circuit_model::bool_to_int(inputs, this->k)];
}






