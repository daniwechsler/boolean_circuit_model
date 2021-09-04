/*
 * circuitmodel.cpp
 */

#include "circuitmodel.h"
#include "booleanfunction.h"
#include <assert.h>
#include <random>
#include <chrono>
#include <vector>
#include <math.h>


circuit_model::circuit_model(vector<boolean_circuit_genome*> active, vector<bool*> passive,
		vector<int> N_active, vector<int> N_passive, boolean_function* f_passive, double b,
		double s, double d, double r_a, double r_p, double alpha, double sigma_a, double sigma_p, double mu_a, double mu_p,
		bool active_sexual_reproduction, bool passive_sexual_reproduction, unsigned long seed, bool emulate_python) {

	assert(this->N_active.size() == this->active.size());
	assert(this->N_passive.size() == this->passive.size());

	this->active = active;
	this->passive = passive;
	this->N_active = N_active;
	this->N_passive = N_passive;
	this->b = b;
	this->s = s;
	this->d = d;
	this->r_a = r_a;
	this->r_p = r_p;
	this->alpha = alpha;
	this->sigma_a = sigma_a;
	this->sigma_p = sigma_p;
	this->mu_a = mu_a;
	this->mu_p = mu_p;
	this->active_sexual_reproduction = active_sexual_reproduction;
	this->passive_sexual_reproduction = passive_sexual_reproduction;
	this->f_passive = f_passive;
	this->seed = seed;
	this->num_traits = this->active[0]->n-2;

	assert(this->alpha <= 1.0 and this->alpha > 0.0);
	assert(this->sigma_a < 1.0 and this->sigma_a >= 0.0);
	assert(this->sigma_p < 1.0 and this->sigma_p >= 0.0);
	assert(this->mu_a < 1.0 and this->mu_a >= 0.0);
	assert(this->mu_p < 1.0 and this->mu_p >= 0.0);

	this->N_a = 0;
	for (unsigned int i=0; i<this->N_active.size(); i++) {
		this->N_a += this->N_active.at(i);
	}
	this->N_p = 0;
	for (unsigned int j=0; j<this->N_passive.size(); j++) {
		this->N_p += this->N_passive.at(j);
	}

	this->payoff_active = new double[this->N_active.size()]();
	this->payoff_passive = new double[this->N_passive.size()]();

	this->t = 0;

	if (this->seed == 0) {
		this->seed = chrono::system_clock::now().time_since_epoch().count();
	}
	this->rand_generator = new rand_man(this->seed, emulate_python);
	this->payoff();
}

circuit_model::~circuit_model() {

}


void circuit_model::payoff(void) {

	delete[] this->payoff_active;
	delete[] this->payoff_passive;
	this->payoff_active = new double[this->N_active.size()]();
	this->payoff_passive = new double[this->N_passive.size()]();

	double f_i, s_j;
	double n_j, n_i_N_p;

	// For each active type the fraction of passive individuals it interacts with
	double* c_active = new double[this->N_active.size()]();

	for (unsigned int i=0; i<this->N_active.size(); i++) {

		for (unsigned int j=0; j<this->N_passive.size(); j++) {

			f_i = (double)((int)this->compute(this->active.at(i), this->passive.at(j)));
			n_j = (double)this->N_passive[j] / (double)this->N_p;
			s_j = (double)((int)this->f_passive->compute(this->passive.at(j)));
			this->payoff_active[i] += n_j * f_i * (s_j*this->r_a + (1.0-s_j) * this->s);
			c_active[i] += n_j * f_i;
		}

		this->payoff_active[i] /= (this->alpha + (1.0-this->alpha) * c_active[i]);
		n_i_N_p = (double)this->N_active[i] / (double)this->N_p;

		for (unsigned int j=0; j<this->N_passive.size(); j++) {
			f_i = (double)((int)this->compute(this->active.at(i), this->passive.at(j)));
			s_j = (double)((int)this->f_passive->compute(this->passive.at(j)));
			this->payoff_passive[j] += n_i_N_p * f_i * (s_j * this->r_p + (1.0 - s_j) * this->b) * (1.0 / (this->alpha + (1.0-this->alpha) * c_active[i]));
		}

	}

	for (unsigned int i=0; i<this->N_active.size(); i++) {
		this->payoff_active[i] -= (double)this->computational_costs(i) * this->d;
	}
	delete[] c_active;
}

bool circuit_model::compute(boolean_circuit_genome* active, bool* passive) {

	bool* inputs = new bool[this->num_traits+2]();
	inputs[0] = false;
	inputs[1] = true;
	for (unsigned int g=0; g<this->num_traits; g++) {
		inputs[g+2] = passive[g];
	}
	bool* res = active->compute(inputs);
	bool res0 = res[0];
	delete[] inputs;
	delete[] res;
	return res0;
}

int circuit_model::computational_costs(int i) {

	int cost = 0;
	for (int k=2; k<this->active[i]->n+this->active[i]->r*this->active[i]->c; k++) {
		cost += (int)this->active[i]->incl_nodes[k];
	}
	return cost;
}


void circuit_model::resample_active(void) {

	int max_comp_cost = this->active[0]->r * this->active[0]->c + this->active[0]->n - 2;
	double min_payoff = this->s - this->d * (double)(max_comp_cost);

	double w_tot = 0.0;
	double* w_active = new double[this->N_active.size()]();
	vector<double> p_sel_active;

	for (unsigned int i=0; i<this->N_active.size(); i++) {
		w_active[i] = (1.0 - this->sigma_a) + this->sigma_a * (this->payoff_active[i] - min_payoff);
		w_tot += w_active[i] * (double)this->N_active[i];
	}

	for (unsigned int i=0; i<this->N_active.size(); i++) {
		p_sel_active.push_back(w_active[i] * (double)this->N_active[i] / w_tot);
	}

	int num_parents = this->N_a;
	// If reproduction is sexual there are twice as many parents needed.
	if (this->active_sexual_reproduction) {
		num_parents = this->N_a*2;
	}
	// Select random parents (proportional to fitness)
	vector<int> parents = this->rand_generator->rand_discrete(p_sel_active, num_parents);

	vector<boolean_circuit_genome*> active_new;
	vector<int> N_active_new;

	boolean_circuit_genome* parent;
	boolean_circuit_genome* parent_two;
	boolean_circuit_genome* offspring;

	for (int i=0; i<this->N_a; i++) {

		if (this->active_sexual_reproduction) {
			parent = this->active[parents[i]];
			parent_two = this->active[parents[this->N_a+i]];
			offspring = parent->recombine(parent_two);
		} else {
			parent = this->active[parents[i]];
			offspring = parent->clone();
		}

		offspring = this->mutate_active_offspring(offspring);
		bool on_list = false;
		for (unsigned int k=0; k<active_new.size(); k++) {

			if (active_new[k]->is_equal(offspring)) {
				N_active_new[k] += 1;
				on_list = true;
				delete offspring;
				break;
			}
		}

		if (!on_list) {
			N_active_new.push_back(1);
			active_new.push_back(offspring);
		}
	}

	for(unsigned int i = 0 ; i < this->active.size(); i++) {
	 	delete this->active[i];
	}
	this->active.clear();
	this->N_active.clear();

	this->active.swap(active_new);
	this->N_active.swap(N_active_new);

	delete[] w_active;
}

boolean_circuit_genome* circuit_model::mutate_active_offspring(boolean_circuit_genome* offspring) {

	for (int g=0; g<offspring->length(); g++) {

		if (this->rand_generator->rand_uniform() <= this->mu_a) {
			// Use manual random generator to find new gene value in
			// case python random generator is emulated.
			// Only for testing purposes (comparing python and C++ version of the model)
			if (this->rand_generator->emulate_python) {
				vector<int> allowed;
				allowed = offspring->get_allowed_values(g);
				vector<double> p;
				for (unsigned int v=0; v<allowed.size(); v++) {
					p.push_back(1.0 / (double)allowed.size());
				}
				int value = this->rand_generator->rand_discrete(p, 1)[0];
				offspring->mutate(g, allowed[value]);
			} else {
				offspring->mutate(g, -1);
			}
		}
	}
	return offspring;
}

void circuit_model::resample_passive(void) {

	double min_payoff = 0.0;
	double w_tot = 0.0;
	double* w_passive = new double[this->N_passive.size()]();
	vector<double> p_sel_passive;

	// Compute absolute fitness
	for (unsigned int j=0; j<this->N_passive.size(); j++) {
		w_passive[j] = (1.0 - this->sigma_p) + this->sigma_p * (this->payoff_passive[j] - min_payoff);
		w_tot += w_passive[j] * (double)this->N_passive[j];
	}

	// Compute selection probability
	for (unsigned int j=0; j<this->N_passive.size(); j++) {
		p_sel_passive.push_back(w_passive[j] * (double)this->N_passive[j] / w_tot);
	}

	vector<bool*> passive_new;
	vector<int> N_passive_new;

	int num_parents = this->N_p;
	// If reproduction is sexual there are twice as many parents needed.
	if (this->passive_sexual_reproduction) {
		num_parents = this->N_p*2;
	}

	vector<int> parents = this->rand_generator->rand_discrete(p_sel_passive, num_parents);

	bool* offspring;
	bool* parent;
	bool* parent_two;

	for (int j=0; j<this->N_p; j++) {

		if (this->passive_sexual_reproduction) {
			parent = this->passive[parents[j]];
			parent_two = this->passive[parents[this->N_p+j]];
			offspring = this->recombine_passive(parent, parent_two);
		} else {
			parent = this->passive[parents[j]];
			offspring = this->clone_passive(parent);
		}
		offspring = this->mutate_passive_offspring(offspring);

		bool on_list = false;
		for (unsigned int k=0; k<passive_new.size(); k++) {
			bool eq = true;
			for (unsigned int g=0; g<this->num_traits; g++) {
				if (offspring[g] != passive_new[k][g]) {
					eq = false;
				}
			}
			if (eq) {
				N_passive_new[k] += 1;
				on_list = true;
				delete[] offspring;
				break;
			}
		}
		if (!on_list) {
			N_passive_new.push_back(1);
			passive_new.push_back(offspring);
		}
	}

	for(unsigned int j = 0 ; j < this->passive.size(); j++) {
	   delete[] this->passive[j];
	}
	this->passive.clear();
	this->N_passive.clear();

	this->passive.swap(passive_new);
	this->N_passive.swap(N_passive_new);
	delete[] w_passive;
}

bool* circuit_model::clone_passive(bool* parent) {
	bool* clone = new bool[this->num_traits]();
	for (unsigned int g=0; g<this->num_traits; g++) {
		clone[g] = parent[g];
	}
	return clone;
}

bool* circuit_model::recombine_passive(bool* parent_one, bool* parent_two) {

	bool* recombinant = new bool[this->num_traits]();
	int crossing_over_point = rand() % this->num_traits;

	for (unsigned int g=0; g<this->num_traits; g++) {
		if (g < crossing_over_point) {
			recombinant[g] = parent_one[g];
		} else {
			recombinant[g] = parent_two[g];
		}
	}
	return recombinant;
}

bool* circuit_model::mutate_passive_offspring(bool* offspring) {

	for (unsigned int g=0; g<this->num_traits; g++) {

		if (this->rand_generator->rand_uniform() <= this->mu_p) {
			if (this->rand_generator->rand_uniform() <= 0.5) {
				offspring[g] = true;
			} else {
				offspring[g] = false;
			}
		}
	}
	return offspring;
}


void circuit_model::sample_time_point(int t, int sample_index) {

	this->T.push_back(t);

	// Sample passive individual data
	vector<int> N_passive_t;
	vector<double> payoff_passive_t;
	for (int j=0; j < pow(2, this->num_traits); j++) {
		N_passive_t.push_back(0);
		payoff_passive_t.push_back(0);
	}

	for (unsigned int j=0; j < this->passive.size(); j++) {
		int dec = circuit_model::bool_to_int(this->passive[j], this->num_traits);
		N_passive_t[dec] += this->N_passive[j];
		payoff_passive_t[dec] += this->payoff_passive[j];
	}

	this->N_passive_T.push_back(N_passive_t);
	this->payoff_passive_T.push_back(payoff_passive_t);

	// Sample active individual data
	unordered_map<unsigned long long, vector<int>>::iterator cursor;
	for (cursor = this->N_active_T.begin(); cursor != this->N_active_T.end(); cursor++) {
		this->N_active_T[cursor->first].push_back(0);
		this->avg_gates_T[cursor->first].push_back(0.0);
		this->payoff_active_T[cursor->first].push_back(0.0);
	}

	for (unsigned int i=0; i < this->active.size(); i++) {

		bool* pehotype_bin = circuit_model::compute_phenotype(this->active[i]);

		unsigned long long phenotype_dec = circuit_model::bool_to_int(pehotype_bin, pow(2, this->num_traits));
		delete[] pehotype_bin;
		double comp_cost = (double)computational_costs(i);

		if (this->N_active_T.find(phenotype_dec) == this->N_active_T.end()) {
			vector<int> v;
			vector<double> w;
			vector<double> x;
			this->N_active_T[phenotype_dec] = v;
			this->avg_gates_T[phenotype_dec] = w;
			this->payoff_active_T[phenotype_dec] = x;
		}
		for (int k=this->N_active_T[phenotype_dec].size(); k<sample_index+1; k++) {
			this->N_active_T[phenotype_dec].push_back(0);
			this->avg_gates_T[phenotype_dec].push_back(0.0);
			this->payoff_active_T[phenotype_dec].push_back(0.0);
		}
		this->N_active_T[phenotype_dec][sample_index] += this->N_active[i];
		this->avg_gates_T[phenotype_dec][sample_index] += (comp_cost * (double)this->N_active[i]);
		this->payoff_active_T[phenotype_dec][sample_index] += (this->payoff_active[i] * (double)this->N_active[i]);
	}

	for (cursor = this->N_active_T.begin(); cursor != this->N_active_T.end(); cursor++) {
		this->avg_gates_T[cursor->first][sample_index] /= (double)this->N_active_T[cursor->first][sample_index];
		this->payoff_active_T[cursor->first][sample_index] /= (double)this->N_active_T[cursor->first][sample_index];
	}
}


void circuit_model::update(void) {

	this->resample_active();
	this->resample_passive();
	this->t += 1;
	this->payoff();
}


void circuit_model::run(int T, int SAMPLE_INTERVAL, int SAMPLE_START, int SAMPLE_END) {

	int sample_index = 0;
	for (int t=0; t<T; t++) {

		if (t >= SAMPLE_START && t<=SAMPLE_END &&  (t % SAMPLE_INTERVAL)==0) {

			this->sample_time_point(t, sample_index);
			sample_index += 1;
		}
		this->update();
	}
}


void circuit_model::run_increment_b(int T, int SAMPLE_INTERVAL, int SAMPLE_START, int SAMPLE_END, int b_level_time, double delta_b, int burn_in) {

	// Update for burn_in number of time steps. In the results table these iterations
	// are not counted as increments of t.
	for (int t=0; t<burn_in; t++) {
		this->update();
	}

	int sample_index = 0;
	for (int t=0; t<T; t++) {

		if(t>0 && t%b_level_time == 0) {
			this->b += delta_b;
		}

		if (t >= SAMPLE_START && t<=SAMPLE_END &&  (t % SAMPLE_INTERVAL)==0) {
			this->sample_time_point(t, sample_index);
			sample_index += 1;
		}
		this->update();
	}
}


bool* circuit_model::compute_phenotype(boolean_circuit_genome* genome) {

	bool* pehnotype = new bool[(int)pow(2, this->num_traits)]();
	for (int j=0; j<pow(2, this->num_traits); j++) {
		bool* inp = circuit_model::int_to_bool(j, this->num_traits);
		pehnotype[j] = this->compute(genome, inp);
		delete[] inp;
	}
	return pehnotype;
}


void circuit_model::print_passive(bool* passive) {
	for (unsigned int g=0; g<this->num_traits; g++) {
		cout << passive[g];
	}
	cout << "\n";
}



