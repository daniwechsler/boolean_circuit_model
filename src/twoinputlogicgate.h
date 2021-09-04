/*
 * logicgate.h
 */

#ifndef TWOINPUTLOGICGATE_H_
#define TWOINPUTLOGICGATE_H_

typedef bool (*fp)(bool, bool);

class two_input_logic_gate {

public:
	two_input_logic_gate(fp func);
	virtual ~two_input_logic_gate();

	int arity;
	fp func;
	bool compute(bool* inputs);

	static bool logic_and(bool x, bool y) {
		return x && y;
	}

	static bool logic_or(bool x, bool y) {
		return x || y;
	}

	static bool logic_nor(bool x, bool y) {
		return !(x || y);
	}

	static bool logic_nand(bool x, bool y) {
		return !(x && y);
	}

	static bool logic_xor(bool x, bool y) {
		return (x || y) && (!(x && y));
	}

	static bool logic_xnor(bool x, bool y) {
		return !((x || y) && (!(x && y)));
	}

	static bool logic_not(bool x, bool y) {
		return !x;
	}
};

#endif /* TWOINPUTLOGICGATE_H_ */
