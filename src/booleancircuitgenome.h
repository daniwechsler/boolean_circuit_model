/*
 * booleancircuitgenome.h
 */

#ifndef BOOLEANCIRCUITGENOME_H_
#define BOOLEANCIRCUITGENOME_H_

#include <iostream>
#include <vector>
#include "twoinputlogicgate.h"

using namespace std;

class boolean_circuit_genome {

public:

	boolean_circuit_genome(int* G, int n, int m, int r, int c, int l, std::vector<two_input_logic_gate*> F);
	virtual ~boolean_circuit_genome();

	int gene(int j);
	int col(int i);
	void decode(void);
	void recursive_include(int j);
	bool* compute(bool* inputs);
	void mutate(int i, int value);
	vector<int> get_allowed_values(int i);
	boolean_circuit_genome* clone(void);
	boolean_circuit_genome* recombine(boolean_circuit_genome* other_genome);
	boolean_circuit_genome* get_random_genome(void);
	int length(void);
	void assert_valid(void);
	bool is_equal(boolean_circuit_genome* other_genome);

	int n;
	int m;
	int c;
	int r;
	int l;
	int* G;
	int arity;
	std::vector<two_input_logic_gate*> F;
	bool* incl_nodes;

	/**
	 * Returns a genome that contains only 0's.
	 */
	static int* create_blank_genome(int r, int c, int m, int arity) {
		int length = (c*r)*(arity+1)+m;
		int* G = new int[length]();
		for (int i=0; i<length; i++) {
			G[i] = 0;
		}
		return G;
	}

};

#endif /* BOOLEANCIRCUITGENOME_H_ */
