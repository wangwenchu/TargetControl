#include "Graph.h"
#include "Util.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <assert.h>

using namespace std;


int main() {
	srand(time(NULL));
	bool isTargetControl = true;
	//isTargetControl = false;
	string gpath = "ieee2224_2.txt";
	string outfile = "driver_node_frequcency.txt";
	Graph link;
	link.read_edgelist(gpath);

	if (isTargetControl) {
		bool isTargetSetGiven = true;
		if (isTargetSetGiven)
		{
			string target_set_path = "fu.txt";
			link.read_target_set(target_set_path);
		}
		else {
			double fraction = 0.2;
			link.random_select_target_node(fraction);
		}

		int sampling_times = 1000 * link.get_target_set_size();
		vector<int>matching_count = link.target_control_node_sampling(sampling_times);
		compute_frequence_and_write(outfile, matching_count, sampling_times);
	}
	else {
		int sampling_times = 1000 * link.get_num_of_nodes();
		vector<int>matching_count = link.full_contrl_node_sampling(sampling_times);
		compute_frequence_and_write(outfile, matching_count, sampling_times);
	}
	std::system("pause");
}