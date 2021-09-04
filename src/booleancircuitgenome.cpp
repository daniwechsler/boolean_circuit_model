/*
 * booleancircuitgenome.cpp
 */

#include "booleancircuitgenome.h"
#include "twoinputlogicgate.h"
#include <assert.h>
#include <stdlib.h>

boolean_circuit_genome::boolean_circuit_genome(int* G, int n, int m, int r, int c, int l, std::vector<two_input_logic_gate*> F) {

	this->G = G;
	this->n = n;
	this->m = m;
	this->r = r;
	this->c = c;
	this->l = l;
	this->F = F;
	this->arity = F.at(0)->arity;
	this->incl_nodes = new bool[this->n+this->r*this->c+this->m]();

	this->assert_valid();
	this->decode();
}

boolean_circuit_genome::~boolean_circuit_genome() {
	delete[] this->G;
	delete[] this->incl_nodes;
}


int boolean_circuit_genome::length(void) {
	return (this->r*this->c)*(this->arity+1) + this->m;
}


int boolean_circuit_genome::gene(int j) {

	return (j-this->n) * (1+this->arity);
}


int boolean_circuit_genome::col(int i) {

	if (this->r == 0) {
		assert(i < this->m);
		return 0;
	}

	return (i-(i % (this->r*(1+this->arity)))) / ((1+this->arity)*this->r);
}


void boolean_circuit_genome::decode(void) {

	for (int j=0; j<this->n+this->r*this->c; j++) {
		this->incl_nodes[j] = false;
	}
	for (int j=this->n+this->r*this->c; j<this->n+this->r*this->c+this->m; j++) {
		this->incl_nodes[j] = true;
	}

	// Start recursive inclusion from reach output node.
	for (int i=this->length()-this->m; i<this->length(); i++) {
		this->recursive_include(this->G[i]);
	}
}


void boolean_circuit_genome::recursive_include(int j) {

	if (j < this->n) {
		this->incl_nodes[j] = true;
		return;
	}
	if (this->incl_nodes[j]) return;

	this->incl_nodes[j] = true;
	for (int i=this->gene(j)+1; i<this->gene(j)+this->arity+1; i++) {
		this->recursive_include(this->G[i]);
	}
}


bool* boolean_circuit_genome::compute(bool* inputs) {

	// Array holds for each gate its output
	bool* outputs = new bool[this->n+this->r*this->c+this->m]();

	// Set outputs of input nodes (just copy)
	for (int i=0; i<this->n; i++) {
		outputs[i] = inputs[i];
	}

	// Compute outputs of internal nodes
	bool* gate_inputs = new bool[this->arity]();
	int gate_input_index = 0;
	two_input_logic_gate* f;
	for (int j=this->n; j<this->n+this->c*this->r; j++) {
		if (this->incl_nodes[j]) {

			f = this->F.at(this->G[this->gene(j)]);

			gate_input_index = 0;
			for (int i=this->gene(j)+1; i<this->gene(j)+this->arity+1; i++) {
				gate_inputs[gate_input_index] = outputs[this->G[i]];
				gate_input_index += 1;
			}
			outputs[j] = f->compute(gate_inputs);
		}
	}
	delete[] gate_inputs;
	// Compute outputs of output nodes
	int j;
	int i;
	bool* external_outputs = new bool[this->m]();
	for (int k=0; k<this->m; k++) {
		i = this->r*this->c*(1+this->arity) + k;
		j = this->n + this->r*this->c + k;
		outputs[j] = outputs[this->G[i]];
		external_outputs[k] = outputs[j];
	}
	delete[] outputs;
	return external_outputs;
}


vector<int> boolean_circuit_genome::get_allowed_values(int i) {

	vector<int> allowed_values;
	if (i % (1 + this->arity) == 0 && i < this->length()-this->m) {
		for (int j=0; j<this->F.size(); j++) {
			allowed_values.push_back(j);
		}
	} else {
		if (this->col(i) >= this->l) {
			// External inputs are allowed
			for (int j=0; j<this->n; j++) {
				allowed_values.push_back(j);
			}
			// Internal gates
			for (int j=this->col(i) * this->r + this->n - this->l * this->r; j<this->col(i) * this->r + this->n; j++) {
				allowed_values.push_back(j);
			}
		} else {
			// All previous gates are allowed
			for (int j=0; j<this->col(i)*this->r+this->n; j++) {
				allowed_values.push_back(j);
			}
		}
	}
	return allowed_values;
}


void boolean_circuit_genome::assert_valid(void) {

	for (unsigned int j=0; j<this->F.size(); j++) {
		assert(this->F.at(j)->arity == this->arity);
	}

	assert(this->length() == (this->r*this->c*(this->arity+1) + this->m));

	vector<int> allowed;
	bool found;
	for (int i=0; i<this->length(); i++) {
		allowed = this->get_allowed_values(i);
		found = false;
		for (unsigned int k=0; k<allowed.size(); k++) {
			if (this->G[i] == allowed.at(k)) found = true;
		}
		assert(found);
	}
}


boolean_circuit_genome* boolean_circuit_genome::clone(void) {

	int* G = new int[this->length()]();
	for (int i=0; i<this->length(); i++) {
		G[i] = this->G[i];
	}
	boolean_circuit_genome* clone = new boolean_circuit_genome(G, this->n, this->m, this->r, this->c, this->l, this->F);
	return clone;
}


boolean_circuit_genome* boolean_circuit_genome::recombine(boolean_circuit_genome* other_genome) {
	/**
	 * Returns a new boolean_circuit_genome which is a recombination of this genome
	 * and the genome given as a parameter.
	 *
	 * The method selects a single crossing over point x at random (an index of the G array).
	 * A position i of G in the recombinant genome contains the value this->G[i] if
	 * i < x otherwise the value other_genome->G[i]
	 */

	assert(other_genome->length() == this->length());

	// Chose random crossing over point
	int crossing_over_point = rand() % this->length();

	// Do crossing over
	int* G = new int[this->length()]();

	for (int i=0; i<this->length(); i++) {
		if (i < crossing_over_point) {
			G[i] = this->G[i];
		} else {
			G[i] = other_genome->G[i];
		}
	}

	boolean_circuit_genome* recombinant = new boolean_circuit_genome(G, this->n, this->m, this->r, this->c, this->l, this->F);
	return recombinant;
}


void boolean_circuit_genome::mutate(int i, int value) {

	assert(i < this->length());
	vector<int> allowed;
	allowed = this->get_allowed_values(i);
	if (value != -1) {
		bool value_valid = false;
		for (unsigned int v=0; v<allowed.size(); v++) {
			if (allowed[v] == value) value_valid = true;
		}
		assert(value_valid);
		this->G[i] = value;
	} else {
		this->G[i] = allowed.at(rand() % allowed.size());
	}
	this->decode();
}


boolean_circuit_genome* boolean_circuit_genome::get_random_genome(void) {

	int* G = boolean_circuit_genome::create_blank_genome(this->r, this->c, this->m, this->arity);
	vector<int> allowed;

	for (int i=0; i<this->length(); i++) {
		allowed = this->get_allowed_values(i);
		G[i] = allowed.at(rand() % allowed.size());
	}
	boolean_circuit_genome* random = new boolean_circuit_genome(G, this->n, this->m, this->r, this->c, this->l, this->F);
	return random;
}


bool boolean_circuit_genome::is_equal(boolean_circuit_genome* other_genome) {

	if (this->length() != other_genome->length()) return false;
	for (int g=0; g<this->length(); g++) {
		if (other_genome->G[g] != this->G[g]) return false;
	}
	return true;
}


