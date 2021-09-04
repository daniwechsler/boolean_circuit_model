/*
 * rand.h
 */

#ifndef RAND_H_
#define RAND_H_

#include <gsl/gsl_rng.h>
#include <vector>

using namespace std;

class rand_man
{
public:

	rand_man(unsigned long int seed, bool emulate_python);
	double rand_uniform();
	vector<int> rand_discrete(vector<double> p, unsigned int size);
	virtual ~rand_man();

	bool emulate_python;
	unsigned long int seed;
	gsl_rng* rng;
};

#endif /* RAND_H_ */


