/*
 * rand.cpp
 */

#include "randman.h"
#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>

using namespace std;

rand_man::rand_man(unsigned long int seed, bool emulate_python) {
	this->emulate_python = emulate_python;
	this->seed = seed;
	this->rng = gsl_rng_alloc( gsl_rng_mt19937 );
	gsl_rng_set(this->rng, this->seed);
}

rand_man::~rand_man(){}

double rand_man::rand_uniform() {
	double r = gsl_rng_uniform( this->rng );
	// If this switch is on, two random numbers are drawn.
	// The second one is not used.
	// If the switch is on repeated calls to rand_uniform will produce
	// the same sequence as repeated calls to np.random.RandomState(seed).random_sample()
	// in python and thus allows to compare c++ with python implementation.
	if (this->emulate_python) {
		gsl_rng_uniform( this->rng );
	}
	return r;
}

vector<int> rand_man::rand_discrete(vector<double> p, unsigned int size) {

	vector<double> F;
	F.push_back(0.0);
	vector<int> samples;
	double sum = 0.0;

	for (unsigned int i=0; i<p.size(); i++) {
		F.push_back(0);
		F[i+1] = F[i] + p[i];
		sum += p[i];
	}
	if (abs(sum-1.0) > 1e-8) {
		assert(false);
	}
	for (unsigned int s=0; s<size; s++) {
		double r = this->rand_uniform();
		for (unsigned int i=0; i<p.size(); i++) {
			if (F[i] <= r and r < F[i+1]) {
				samples.push_back(i);
				break;
			}
		}
	}
	return samples;
}
