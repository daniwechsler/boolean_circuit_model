/*
 * circuitmodel.h
 */

#ifndef SRC_CIRCUITMODEL_H_
#define SRC_CIRCUITMODEL_H_

#include "booleancircuitgenome.h"
#include "booleanfunction.h"
#include "randman.h"
#include <vector>
#include <random>
#include <unordered_map>

/**
 * A passive individual is implemented as an array of booleans.
 *
 * position: [0  1  2]
 * value:	 [T  T  F]
 *
 *
 * The decimal representation of a passive individual is obtained
 * by converting the binary representation to a decimal number with
 * the low position bits being the most significant bits.
 *
 * bin([T T F]) = 1*2^2 + 1*2^1 + 0*2^0 6
 *
 *
 */


class circuit_model {

public:

	circuit_model(vector<boolean_circuit_genome*> active, vector<bool*> passive,
			vector<int> N_active, vector<int> N_passive, boolean_function* f_passive, double b, double s, double d,
			double r_a, double r_p, double alpha, double sigma_a, double sigma_p, double mu_a, double mu_p,
			bool active_sexual_reproduction, bool passive_sexual_reproduction, unsigned long seed, bool emulate_python);
	virtual ~circuit_model();

	void payoff(void);
	int computational_costs(int i);
	bool compute(boolean_circuit_genome* active, bool* passive);
	void resample_active(void);
	boolean_circuit_genome* mutate_active_offspring(boolean_circuit_genome* parent);
	void resample_passive(void);
	bool* clone_passive(bool* parent);
	bool* recombine_passive(bool* parent_one, bool* parent_two);
	bool* mutate_passive_offspring(bool* parent);
	bool* compute_phenotype(boolean_circuit_genome* genome);
	void sample_time_point(int t, int sample_index);
	void update(void);
	void run(int T, int SAMPLE_INTERVAL, int SAMPLE_START, int SAMPLE_END);
	void run_increment_b(int T, int SAMPLE_INTERVAL, int SAMPLE_START, int SAMPLE_END, int b_level_time, double delta_b, int burn_in);
	void print_passive(bool* passive);

	vector<boolean_circuit_genome*> active;
	vector<bool*> passive;
	vector<int> N_active;
	vector<int> N_passive;
	boolean_function* f_passive;
	double b;
	double s;
	double d;
	double r_a;
	double r_p;
	double alpha;
	double sigma_a;
	double sigma_p;
	double mu_a;
	double mu_p;
	bool active_sexual_reproduction;
	bool passive_sexual_reproduction;
	double* payoff_active;
	double* payoff_passive;
	int N_a;
	int N_p;
	int t;
	unsigned int num_traits;
	unsigned int seed;
	rand_man* rand_generator;

	vector<int> T;
	vector<vector<int>> N_passive_T;
	vector<vector<double>> payoff_passive_T;
	unordered_map<unsigned long long, vector<int>> N_active_T;
	unordered_map<unsigned long long, vector<double>> avg_gates_T;
	unordered_map<unsigned long long, vector<double>> payoff_active_T;

	/**
	 * Converts an array of boolean to a decimal representation.
	 * The bit with the smallest index is the most significant bit.
	 *
	 * Example:
	 * bool_to_int([0 1 1 1], 4) = 7
	 */
	static unsigned long bool_to_int(bool* b, unsigned int size) {
		unsigned long long sum = 0;
		for (unsigned int i=0; i<size; i++) {
			sum += long(b[i]) * pow(2, size-i-1);
		}
		return sum;
	}

	/**
	 * Converts a given integer 'num' to its binary representation.
	 * The binary representation is an array of boolean values with
	 * the values with low indices being the most significant
	 * bits.
	 *
	 * Example:
	 * int_to_bool(7, 4) = [0 1 1 1]
	 */
	static bool* int_to_bool(unsigned int num, unsigned int size) {
		bool* b = new bool[size]();

		for (unsigned int i = 0; i < size; i++)
		{
		    b[size-i-1] = ((num >> i) & 1) == 1;
		}
		return b;
	}


};

#endif /* SRC_CIRCUITMODEL_H_ */
