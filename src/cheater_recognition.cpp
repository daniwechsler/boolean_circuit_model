#include <iostream>
#include "booleancircuitgenome.h"
#include "twoinputlogicgate.h"
#include "booleanfunction.h"
#include "circuitmodel.h"
#include "randman.h"
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <tclap/CmdLine.h>
#include <string>
#include <time.h>

using namespace std;

int T;
int SAMPLE_INTERVAL;
int SAMPLE_START;
int SAMPLE_END;
int RANDOM_SEED;
int n;
int m;
int c;
int r;
int l;
string omega_str;
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
int N_a;
int N_p;
string f_passive_str;

int b_level_time;
double delta_b;
int burn_in;

bool verbose;
bool emulate_python;

void initCmdArguments ( int argc, char* argv[], TCLAP::CmdLine *cmd)
{

	TCLAP::ValueArg<int> T_arg("", "T", "Number of generations to simulate (default is 100).", false, 100, "integer");
	TCLAP::ValueArg<int> SAMPLE_INTERVAL_arg("", "SAMPLE_INTERVAL", "Sample interval (default is 10).", false, 10, "integer");
	TCLAP::ValueArg<int> SAMPLE_START_arg("", "SAMPLE_START", "Start sampling time.", false, 0, "integer");
	TCLAP::ValueArg<int> SAMPLE_END_arg("", "SAMPLE_END", "End sampling time.", false, 0, "integer");
	TCLAP::ValueArg<int> RANDOM_SEED_arg("", "RANDOM_SEED", "Random seed.", false, 0, "integer");

	TCLAP::ValueArg<int> n_arg("", "n", "Number of circuit inputs (default is 5)", false, 5, "integer");
	TCLAP::ValueArg<int> c_arg("", "c", "Number of circuit columns (default is 2).", false, 2, "integer");
	TCLAP::ValueArg<int> r_arg("", "r", "Number of circuit rows (default is 2).", false, 2, "integer");
	TCLAP::ValueArg<int> l_arg("", "l", "Allowed input gates (default is 1)", false, 1, "integer");
	TCLAP::ValueArg<std::string> omega_arg("", "omega","Set of available logic gates (e.g. AND,OR) (default is NAND)", false, "NAND", "string");

	TCLAP::ValueArg<double> b_arg("", "b", "Payoff a cheating plant receives when interaction with an animal (default is 2.0).", false, 2.0, "double");
	TCLAP::ValueArg<double> s_arg("", "s", "Payoff an animal receives when interacting with a cheater (default is -1.0).", false, -1.0, "double");
	TCLAP::ValueArg<double> d_arg("", "d", "Cost an animal pays for each logic gate and input used (default is 0.0).", false, 0.0, "double");
	TCLAP::ValueArg<double> r_a_arg("", "r_a", "Reward the animal gets when both cooperate (default is 1.0).", false, 1.0, "double");
	TCLAP::ValueArg<double> r_p_arg("", "r_p", "Reward the plant gets when both cooperate (default is 1.0).", false, 1.0, "double");
	TCLAP::ValueArg<double> alpha_arg("", "alpha", "Degree to which animals only encounter passive individuals they interact with (default is 1.0).", false, 1.0, "double");

	TCLAP::ValueArg<std::string> f_passive_arg("", "f_passive","Epistatic function given as a sequence of length 2^n consisting of 0's and 1's (default is 11111111).", false, "11111111", "string");

	TCLAP::ValueArg<double> sigma_a_arg("", "sigma_a", "Strength of selection for animals (default is 0.5).", false, 0.5, "double");
	TCLAP::ValueArg<double> sigma_p_arg("", "sigma_p", "Strength of selection for plants (default is 0.5).", false, 0.5, "double");
	TCLAP::ValueArg<double> mu_a_arg("", "mu_a", "Mutation rate of animals (default is 0.01).", false, 0.01, "double");
	TCLAP::ValueArg<double> mu_p_arg("", "mu_p", "Mutation rate of plants (default is 0.01).", false, 0.01, "double");
	TCLAP::SwitchArg active_sexual_reproduction_switch("", "animals_sexual", "If set animals reproduce sexually (default is asexual reproduction).", *cmd, false);
	TCLAP::SwitchArg passive_sexual_reproduction_switch("", "plants_sexual", "If set plants reproduce sexually (default is asexual reproduction).", *cmd, false);

	TCLAP::ValueArg<int> N_a_arg("", "N_a", "Number of animals (default is 100).", false,100, "integer");
	TCLAP::ValueArg<int> N_p_arg("", "N_p", "Number of plants (default is 100).", false,100, "integer");

	TCLAP::SwitchArg verboseSwitch("v", "verbose", "Print the status at every generation.", *cmd, false);
	TCLAP::SwitchArg emulate_python_switch("", "emulate_python", "Whether python random generator should be emulated.", *cmd, false);

	TCLAP::ValueArg<int> b_level_time_arg("", "b_level_time", "b level time.", false, 1000, "integer");
	TCLAP::ValueArg<double> delta_b_arg("", "delta_b", "Delta b.", false, 0.0, "double");
	TCLAP::ValueArg<int> burn_in_arg("", "burn_in", "burn in time.", false, 0, "integer");


	cmd->add(T_arg);
	cmd->add(SAMPLE_INTERVAL_arg);
	cmd->add(SAMPLE_START_arg);
	cmd->add(SAMPLE_END_arg);
	cmd->add(RANDOM_SEED_arg);

	cmd->add(n_arg);
	cmd->add(c_arg);
	cmd->add(r_arg);
	cmd->add(l_arg);
	cmd->add(omega_arg);

	cmd->add(b_arg);
	cmd->add(s_arg);
	cmd->add(d_arg);
	cmd->add(r_a_arg);
	cmd->add(r_p_arg);
	cmd->add(alpha_arg);

	cmd->add(f_passive_arg);

	cmd->add(sigma_a_arg);
	cmd->add(sigma_p_arg);
	cmd->add(mu_a_arg);
	cmd->add(mu_p_arg);

	cmd->add(N_a_arg);
	cmd->add(N_p_arg);

	cmd->add(b_level_time_arg);
	cmd->add(delta_b_arg);
	cmd->add(burn_in_arg);

	cmd->parse( argc, argv );


	T = T_arg.getValue();
	SAMPLE_INTERVAL = SAMPLE_INTERVAL_arg.getValue();
	SAMPLE_START = SAMPLE_START_arg.getValue();
	SAMPLE_END = SAMPLE_END_arg.getValue();
	if (SAMPLE_END == 0) SAMPLE_END = T;
	RANDOM_SEED = RANDOM_SEED_arg.getValue();
	if (RANDOM_SEED == 0) RANDOM_SEED = time(NULL);
	n = n_arg.getValue();
	m = 1;
	c = c_arg.getValue();
	r = r_arg.getValue();
	l = l_arg.getValue();
	omega_str = omega_arg.getValue();
	b = b_arg.getValue();
	s = s_arg.getValue();
	d = d_arg.getValue();
	r_a = r_a_arg.getValue();
	r_p = r_p_arg.getValue();
	alpha = alpha_arg.getValue();
	f_passive_str = f_passive_arg.getValue();
	sigma_a = sigma_a_arg.getValue();
	sigma_p = sigma_p_arg.getValue();
	mu_a = mu_a_arg.getValue();
	mu_p = mu_p_arg.getValue();
	active_sexual_reproduction = active_sexual_reproduction_switch.getValue();
	passive_sexual_reproduction = passive_sexual_reproduction_switch.getValue();
	N_a = N_a_arg.getValue();
	N_p = N_p_arg.getValue();
	b_level_time = b_level_time_arg.getValue();
	delta_b = delta_b_arg.getValue();
	burn_in = burn_in_arg.getValue();

	verbose = verboseSwitch.getValue();
	emulate_python = emulate_python_switch.getValue();

	if (n > 8) {
		// The problem of n>8 concerns only the binary representation of the active phenotype
		// in the output table and not the simulation itself.
		cout << "Larger n leads to overflow of active individual phenotype decimal representation.\n";
		assert(n<=7);
	}

}

void print_rule_table(boolean_circuit_genome* genome) {

	bool* inp;
	bool* outputs;
	bool* inp_r = new bool[genome->n]();
	int k;
	for (int i=0; i<pow(2,genome->n); i++) {
		inp = circuit_model::int_to_bool(i, genome->n);
		k = 0;
		for (int j=genome->n-1; j>=0; j--) {
			cout << inp[j] << " ";
			inp_r[k] = inp[j];
			k += 1;
		}
		outputs = genome->compute(inp_r);
		cout << " : ";
		for (int l=0; l<genome->m; l++) {
			cout << outputs[l] << " ";
		}
		cout << "\n";
	}
}

void print_gnome(boolean_circuit_genome* genome) {

	for (int i=0; i<genome->length(); i++) {
		cout << genome->G[i] << " ";
	}
	cout << "\n";
}

void print_args() {

	cout << "T: \t" << T << "\n";
	cout << "SAMPLE_INTERVAL: \t" << SAMPLE_INTERVAL << "\n";
	cout << "SAMPLE_START: \t" << SAMPLE_START << "\n";
	cout << "SAMPLE_END: \t" << SAMPLE_END << "\n";
	cout << "n: \t" << n << "\n";
	cout << "c: \t" << c << "\n";
	cout << "r: \t" << r << "\n";
	cout << "l: \t" << l << "\n";
	cout << "omega: \t" << omega_str << "\n";
	cout << "b: \t" << b << "\n";
	cout << "s: \t" << s << "\n";
	cout << "d: \t" << d << "\n";
	cout << "r_a: \t" << r_a << "\n";
	cout << "r_p: \t" << r_p << "\n";
	cout << "alpha: \t" << alpha << "\n";
	cout << "sigma_a: \t" << sigma_a << "\n";
	cout << "sigma_p: \t" << sigma_p << "\n";
	cout << "mu_a: \t" << mu_a << "\n";
	cout << "mu_p: \t" << mu_p << "\n";
	cout << "active_sexual_reproduction: \t" << active_sexual_reproduction << "\n";
	cout << "passive_sexual_reproduction: \t" << passive_sexual_reproduction << "\n";
	cout << "N_a: \t" << N_a << "\n";
	cout << "N_p: \t" << N_p << "\n";
	cout << "b_level_time: \t" << b_level_time << "\n";
	cout << "delta_b: \t" << delta_b << "\n";
	cout << "burn_in: \t" << burn_in << "\n";
}

void to_csv(circuit_model* cm, string col_sep) {

	cout << "T" << col_sep;
	for (unsigned int j=0; j<pow(2, cm->num_traits); j++) {
		cout << "N_passive_" << j << col_sep;
	}

	for (unsigned int j=0; j<pow(2, cm->num_traits); j++) {
		cout << "payoff_passive_" << j << col_sep;
	}

	unordered_map<unsigned long long, vector<int>>::iterator cursor;
	vector<unsigned long long> phenotypes;
	for (cursor = cm->N_active_T.begin(); cursor != cm->N_active_T.end(); cursor++) {
		phenotypes.push_back(cursor->first);
	}

	for (unsigned long i=0; i<phenotypes.size(); i++) {
		cout << "N_active_" << phenotypes[i] << col_sep;
	}

	for (unsigned long i=0; i<phenotypes.size(); i++) {
		cout << "avg_payoff_active_" << phenotypes[i] << col_sep;
	}

	for (unsigned long i=0; i<phenotypes.size(); i++) {
		cout << "avg_num_gates_" << phenotypes[i];
		if (i<phenotypes.size()-1) cout << col_sep;
	}

	cout << endl;
	for (unsigned int index=0; index<cm->T.size(); index++) {
		cout << cm->T[index] << col_sep;
		for (unsigned int j=0; j<cm->N_passive_T[index].size(); j++) {
			cout << cm->N_passive_T[index][j] << col_sep;
		}

		for (unsigned int j=0; j<cm->N_passive_T[index].size(); j++) {
			cout << cm->payoff_passive_T[index][j] << col_sep;
		}

		for (unsigned int i=0; i<phenotypes.size(); i++) {
			cout << cm->N_active_T[phenotypes[i]][index]  << col_sep;
		}

		for (unsigned int i=0; i<phenotypes.size(); i++) {
			cout << cm->payoff_active_T[phenotypes[i]][index]  << col_sep;
		}

		for (unsigned int i=0; i<phenotypes.size(); i++) {
			cout << cm->avg_gates_T[phenotypes[i]][index];
			if (i<phenotypes.size()-1) cout << col_sep;
		}

		cout << endl;
	}
}

vector<bool> parse_f_passive_str(string f_passive_str) {

	vector<bool> f_passive;
	for (unsigned int k=0; k<f_passive_str.length(); k++) {
		if (f_passive_str[k] == '1') {
			f_passive.push_back(true);
		} else if (f_passive_str[k] == '0') {
			f_passive.push_back(false);
		} else {
			throw "Invalid f_passive.";
		}
	}
	return f_passive;
}

int main(int argc, char **argv) {


	try {
		TCLAP::CmdLine cmd("Cheater Recognition Model Simulator", ' ', "1.0");
		initCmdArguments(argc, argv, &cmd);

		if (verbose) print_args();

		// Setup initial active agent
		int* G = boolean_circuit_genome::create_blank_genome(r, c, m, 2);

		// Parse list of allowed logic gates
		vector<two_input_logic_gate*> F;

		istringstream iss(omega_str);
		string item;
		while (std::getline(iss, item, ',')) {

			if (item.compare("AND") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_and));
			} else if (item.compare("OR") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_or));
			} else if (item.compare("NOR") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_nor));
			} else if (item.compare("NAND") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_nand));
			} else if (item.compare("XOR") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_xor));
			} else if (item.compare("XNOR") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_xnor));
			} else if (item.compare("NOT") == 0) {
				F.push_back(new two_input_logic_gate(&two_input_logic_gate::logic_not));
			} else {
				cout << "Error: Logic gate '" << item << "' undefined.\n";
				return 1;
			}
		}

		boolean_circuit_genome* active_init = new boolean_circuit_genome(G, n, m, r, c, l, F);

		vector<boolean_circuit_genome*> active;
		vector<int> N_active;
		active.push_back(active_init);
		N_active.push_back(N_a);

		// Setup initial passive agent
		bool* passive_init = new bool[n-2]();
		for (int k=0; k<n-2; k++) {
			passive_init[k] = false;
		}

		vector<bool*> passive;
		vector<int> N_passive;
		passive.push_back(passive_init);
		N_passive.push_back(N_p);

		// Setup f_passive
		vector<bool> f_passive_v = parse_f_passive_str(f_passive_str);
		boolean_function* f_passive = new boolean_function(f_passive_v);

		circuit_model* cm = new circuit_model(active, passive,
				N_active, N_passive, f_passive, b,
				s, d, r_a, r_p, alpha, sigma_a, sigma_p, mu_a, mu_p,
				active_sexual_reproduction, passive_sexual_reproduction, RANDOM_SEED, emulate_python);

		if (delta_b == 0.0) {
			cm->run(T, SAMPLE_INTERVAL, SAMPLE_START, SAMPLE_END);
		} else {
			cm->run_increment_b(T, SAMPLE_INTERVAL, SAMPLE_START, SAMPLE_END, b_level_time, delta_b, burn_in);
		}
		if (verbose) {
			to_csv(cm, "\t");
		} else {
			to_csv(cm, ",");
		}
		return 0;
	}
	catch (TCLAP::ArgException &e) {
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
		return 1;
	}
	catch (std::exception &e) {
		return 1;
	}
}
