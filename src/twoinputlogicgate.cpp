/*
 * logicgate.cpp
 */

#include "twoinputlogicgate.h"

two_input_logic_gate::two_input_logic_gate(fp func) {
	this->arity = 2;
	this->func = func;
}

two_input_logic_gate::~two_input_logic_gate() {
	// TODO Auto-generated destructor stub
}

bool two_input_logic_gate::compute(bool* inputs) {
	return (*this->func)(inputs[0], inputs[1]);
}
