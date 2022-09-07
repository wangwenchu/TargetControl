#include "Graph.h"

Graph::Graph() :num_node(0), num_edge(0)
{
	in_link.clear();
	out_link.clear();
	mark.clear();
	target_set.clear();
	match_left.clear();
	match_right.clear();
}

void Graph::alloc_memory(size_t graphsize)
{
	in_link.resize(graphsize);
	out_link.resize(graphsize);
	mark.resize(graphsize);
	match_left.resize(graphsize);
	match_right.resize(graphsize);
}

bool Graph::left_dfs(int u)
{
	int tmp;
	for (int v : out_link[u]) {
		if (mark[v] == 0) {
			mark[v] = 1;
			tmp = match_right[v];
			if (tmp == -1 || left_dfs(tmp)) {
				match_right[v] = u;
				return true;
			}
		}
	}
	return false;
}

bool Graph::right_dfs(int v)
{
	int tmp;
	for (int u : in_link[v]) {
		if (mark[u] == 0) {
			mark[u] = 1;
			tmp = match_left[u];
			if (tmp == -1 || right_dfs(tmp)) {
				match_left[u] = v;
				return true;
			}
		}
	}
	return false;
}

void Graph::read_edgelist(const string& gpath)
{
	ifstream fin(gpath);
	if (!fin.good()) {
		cerr << "network file error!" << endl;
		std::this_thread::sleep_for(std::chrono::milliseconds(2000)); //单位是毫秒
		exit(1);
	}
	fin >> num_node >> num_edge;
	assert(num_node > 0 && num_edge > 0);
	alloc_memory(num_node);
	int source, target;
	for (int i = 0; i < num_edge; ++i) {
		fin >> source >> target;
		if (source >= num_node || target >= num_node) {
			cerr << "node index error!" << endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(2000)); //单位是毫秒
			exit(1);
		}
		out_link[source].push_back(target);
		in_link[target].push_back(source);
	}
	fin.close();
}

void Graph::read_target_set(const string& target_set_path)
{
	ifstream fin(target_set_path);
	if (!fin.good()) {
		cerr << "target set file error!" << endl;
		std::this_thread::sleep_for(std::chrono::milliseconds(2000)); //单位是毫秒
		exit(1);
	}
	int s;
	while (fin >> s) {
		if (s >= num_node) {
			cerr << "target node error!" << endl;
			std::this_thread::sleep_for(std::chrono::milliseconds(2000)); //单位是毫秒
			exit(1);
		}
		target_set.insert(s);
	}
	fin.close();
}

void Graph::random_select_target_node(double fraction)
{
	int m = num_node * fraction;
	vector<int>node_pool(num_node);
	for (int i = 0; i < num_node; ++i) {
		node_pool[i] = i;
	}
	std::random_shuffle(node_pool.begin(), node_pool.end());
	target_set.clear();
	for (int i = 0; i < m; ++i) {
		target_set.insert(node_pool[i]);
	}
}

std::pair<std::set<int>, std::set<int>> Graph::max_matching(const set<int>& new_target_set)
{
	std::fill(match_left.begin(), match_left.end(), -1);


	set<int>right_match_set;
	for (const int& x : new_target_set) {
		std::fill(mark.begin(), mark.end(), 0);
		if (right_dfs(x)) {
			right_match_set.insert(x);
		}
	}
	set<int>left_match_set;
	std::fill(match_right.begin(), match_right.end(), -1);
	int u;
	for (int i = 0; i < num_node; ++i) {
		u = match_left[i];
		if (u != -1) {
			left_match_set.insert(i);
			match_right[u] = i;
		}
	}

	return { left_match_set, right_match_set };
}

std::pair<std::set<int>, std::set<int>> Graph::full_control()
{
	set<int>cur_target_set;
	for (int i = 0; i < num_node; ++i) {
		cur_target_set.insert(i);
	}
	return max_matching(cur_target_set);
}

std::set<int> Graph::find_left_always_matched_nodes(const vector<int>& cur_match_list, const set<int>& match_set)
{
	//left side
	set<int>always_match_set;
	int node;
	for (int x : match_set) {

		match_left = cur_match_list;
		match_left[x] = -1;
		//cout << x << "\t" << cur_match_list[x] << endl;
		std::fill(mark.begin(), mark.end(), 0);
		mark[x] = 1;
		node = cur_match_list[x];
		if (!right_dfs(node)) {
			always_match_set.insert(x);
		}
	}

	return always_match_set;
}

std::set<int> Graph::find_right_always_matched_nodes(const vector<int>& cur_match_list, const set<int>& match_set, const set<int>& cur_target_set)
{
	set<int>always_match_set;
	int node;
	std::fill(mark.begin(), mark.end(), 1);
	for (int x : match_set) {

		match_right = cur_match_list;
		match_right[x] = -1;
		//cout << x << "\t" << cur_match_list[x] << endl;
		for (int x : cur_target_set) {
			mark[x] = 0;
		}
		mark[x] = 1;
		node = cur_match_list[x];
		if (!left_dfs(node)) {
			always_match_set.insert(x);
		}
	}

	return always_match_set;
}

std::set<int> Graph::get_right_unmatched_nodes(const set<int>& cur_target_set, const set<int>& right_match_set)
{
	set<int>unmatch_set;
	set_difference(cur_target_set.begin(), cur_target_set.end(), right_match_set.begin(), right_match_set.end(), inserter(unmatch_set, unmatch_set.begin()));
	return unmatch_set;
}

std::vector<int> Graph::remove_always_matched_nodes_from_matched_set(const set<int>& matched_set, const set<int>& always_matched)
{
	vector<int>remain_set;
	set_difference(matched_set.begin(), matched_set.end(), always_matched.begin(), always_matched.end(), inserter(remain_set, remain_set.begin()));
	return remain_set;
}

std::vector<int> Graph::find_alternative_node_set(int select_side, int target_node, const vector<int>& cur_match_list)
{
	if (select_side == LEFT) {
		return find_left_altertive_set(target_node, cur_match_list);
	}
	else if (select_side == RIGHT) {
		return find_right_altertive_set(target_node, cur_match_list);
	}
	else {
		cerr << "parameters error!" << endl;
		std::this_thread::sleep_for(std::chrono::milliseconds(2000)); //单位是毫秒
		exit(1);
	}
}

std::vector<int> Graph::find_left_altertive_set(int target_node, const vector<int>& cur_match_list)
{
	/*for (int i = 0; i < cur_match_list.size(); ++i) {
		cout << i << "\t" << cur_match_list[i] << endl;
	}*/
	vector<int>last_match_list;
	vector<int>alternative_set = { target_node };
	int start_node;
	match_left = cur_match_list;  //They might be the same one 
	while (true) {
		last_match_list = match_left;
		start_node = last_match_list[target_node];
		match_left[target_node] = -1;
		std::fill(mark.begin(), mark.end(), 0);
		for (int x : alternative_set) {
			mark[x] = 1;
		}

		if (right_dfs(start_node)) {
			target_node == -1;  //just for test
			for (int i = 0; i < num_node; ++i) {
				if (last_match_list[i] == -1 && match_left[i] != -1) {
					alternative_set.push_back(i);
					target_node = i;
					break;
				}
			}
			assert(target_node != -1);
			last_match_list = match_left;
		}
		else {
			break;
		}
	}
	assert(alternative_set.size() > 1);
	return vector<int>(alternative_set.begin() + 1, alternative_set.end());
}

std::vector<int> Graph::find_right_altertive_set(int target_node, const vector<int>& cur_match_list)
{
	vector<int>last_match_list;
	vector<int>alternative_set = { target_node };
	match_right = cur_match_list; //they might be the same one
	int start_node;
	while (true) {
		last_match_list = match_right;
		start_node = last_match_list[target_node];
		match_right[target_node] = -1;
		std::fill(mark.begin(), mark.end(), 0);
		for (int x : alternative_set) {
			mark[x] = 1;
		}

		if (left_dfs(start_node)) {
			target_node == -1;
			for (int i = 0; i < num_node; ++i) {
				if (last_match_list[i] == -1 && match_right[i] != -1) {
					alternative_set.push_back(i);
					target_node = i;
					break;
				}
			}
			assert(target_node != -1);
			last_match_list = match_left;
		}
		else {
			break;
		}
	}
	assert(alternative_set.size() > 1);
	return vector<int>(alternative_set.begin() + 1, alternative_set.end());
}

std::set<int> Graph::merge_remain_with_always_matched(const vector<int>& remain_set, const set<int>& alway_match)
{
	set<int> match_set = alway_match;
	for (int x : remain_set) {
		match_set.insert(x);
	}
	return match_set;
}

bool Graph::isSubset(const set<int>& small_set, const set<int>& big_set)
{
	if (small_set.empty() || big_set.empty()) return false;
	for (int x : small_set) {
		if (big_set.find(x) == big_set.end()) {
			return false;
		}
	}
	return true;
}

void Graph::remove_repeated_nodes(set<int>& s1, const set<int>& s2)
{
	set<int>::iterator iter;
	for (int x : s2) {
		iter = s1.find(x);
		assert(iter != s1.end());
		s1.erase(iter);
	}
}

void Graph::init_target_node_map_to_original()
{
	target_node_map.clear();
	for (int x : target_set) {
		target_node_map[x] = x;
	}
}

void Graph::update_target_node_map(const set<int>& cur_target_set, const vector<int>& cur_match_right)
{
	int u;
	for (auto& x : target_node_map) {
		u = x.second;
		if (u != -1) {
			x.second = cur_match_right[u];
		}
	}
}

set<int> Graph::update_control_link()
{
	set<int>duplicate_set;
	int u, v;
	for (auto x : target_node_map) {
		u = x.second;
		if (u != -1) {
			v = x.first;
			if (control_link[v].find(u) == control_link[v].end()) {
				control_link[v].insert(u);
			}
			else {
				duplicate_set.insert(u);
			}
		}
	}
	return duplicate_set;
}

std::vector<std::vector<int>> Graph::randomize_node_label_and_link_table(const vector<vector<int>>& link_table)
{
	vector<int>new_node_pool(num_node);    //original(index)-->new(value)
	for (int i = 0; i < num_node; ++i) {
		new_node_pool[i] = i;
	}
	random_shuffle(new_node_pool.begin(), new_node_pool.end());
	vector<int>original_node_pool(num_node);
	for (int i = 0; i < num_node; ++i) {
		original_node_pool[new_node_pool[i]] = i;
	}
	in_link.clear();
	in_link.resize(num_node);
	int node;
	for (int i = 0; i < num_node; ++i) {
		node = new_node_pool[i];
		for (int neibor : link_table[i]) {
			in_link[node].push_back(new_node_pool[neibor]);
		}
		random_shuffle(in_link[node].begin(), in_link[node].end());
	}
	return { new_node_pool,original_node_pool };
}

std::set<int> Graph::get_maped_node(const set<int>& s, const vector<int>& node_pool)
{
	set<int>res;
	for (const int x : s) {
		res.insert(node_pool[x]);
	}
	return res;
}

std::vector<int> Graph::get_maped_node(const vector<int>& s, const vector<int>& node_pool)
{
	vector<int>res(num_node, -1);
	for (int i = 0; i < num_node; ++i) {
		if (s[i] != -1) {
			res[node_pool[i]] = node_pool[s[i]];
		}
	}
	return res;
}

void Graph::unmatch_removed_node_and_max_matching_new_node(int select_side, int unmatch_node, int new_pick_node, vector<int>& cur_match_left, vector<int>& cur_match_right)
{
	if (select_side == LEFT) {
		unmatch_left_node_and_max_matching(unmatch_node, new_pick_node, cur_match_left, cur_match_right);
	}
	else if (select_side == RIGHT) {
		unmatch_right_node_and_max_matching(unmatch_node, new_pick_node, cur_match_left, cur_match_right);
	}
}

void Graph::random_max_matching(set<int>& new_target_set, set<int>& driver_set, const vector<vector<int>>& link_table)
{
	vector<vector<int>> node_pool;
	set<int> new_labeled_target_set;
	pair<set<int>, set<int>> left_right_match_set;
	set<int>left_match_set;
	set<int>right_unmatch_set;
	set<int>duplicate_target_set;
	while (!new_target_set.empty()) {
		node_pool = randomize_node_label_and_link_table(link_table);
		new_labeled_target_set = get_maped_node(new_target_set, node_pool[0]);

		left_right_match_set = max_matching(new_labeled_target_set);
		left_match_set = left_right_match_set.first;
		right_unmatch_set = get_right_unmatched_nodes(new_labeled_target_set, left_right_match_set.second);
		right_unmatch_set = get_maped_node(right_unmatch_set, node_pool[1]);
		driver_set.insert(right_unmatch_set.begin(), right_unmatch_set.end());
		left_match_set = get_maped_node(left_match_set, node_pool[1]);

		match_right = get_maped_node(match_right, node_pool[1]);
		update_target_node_map(new_target_set, match_right);
		duplicate_target_set = update_control_link();
		remove_repeated_nodes(left_match_set, duplicate_target_set);
		new_target_set = left_match_set;
	}
}




void Graph::unmatch_left_node_and_max_matching(int unmatch_node, int new_pick_node, vector<int>& cur_match_left, vector<int>& cur_match_right)
{
	int u = cur_match_left[unmatch_node];

	match_right = cur_match_right;
	match_right[u] = -1;

	bool flag = left_dfs(new_pick_node);
	assert(flag == true && match_right[u] != -1);
	int v;
	cur_match_right = match_right;
	std::fill(cur_match_left.begin(), cur_match_left.end(), -1);
	for (int i = 0; i < num_node; ++i) {
		v = match_right[i];
		if (v != -1) {
			cur_match_left[v] = i;
		}
	}
}

void Graph::unmatch_right_node_and_max_matching(int unmatch_node, int new_pick_node, vector<int>& cur_match_left, vector<int>& cur_match_right)
{
	int u = cur_match_right[unmatch_node];
	match_left = cur_match_left;
	match_left[u] = -1;
	bool flag = right_dfs(new_pick_node);
	assert(flag == true && match_left[u] != -1);
	int v;
	cur_match_left = match_left;
	std::fill(cur_match_right.begin(), cur_match_right.end(), -1);
	for (int i = 0; i < num_node; ++i) {
		v = match_left[i];
		if (v != -1) {
			cur_match_right[v] = i;
		}
	}
}

std::vector<int> Graph::target_control_node_sampling(int sample_times)
{
	printf("%d sampling is beginning\n", sample_times);

	pair<set<int>, set<int>> out_in_match_set = max_matching(target_set);
	set<int> left_match_set = out_in_match_set.first;
	set<int> right_match_set = out_in_match_set.second;

	vector<int> cur_match_left = match_left;
	vector<int> cur_match_right = match_right;
	set<int>left_always = find_left_always_matched_nodes(cur_match_left, left_match_set);
	set<int>right_always = find_right_always_matched_nodes(cur_match_right, right_match_set, target_set);

	vector<int>left_remain = remove_always_matched_nodes_from_matched_set(left_match_set, left_always);
	vector<int>right_remain = remove_always_matched_nodes_from_matched_set(right_match_set, right_always);

	const vector<vector<int>> link_table = in_link;

	vector<int>matching_count(num_node, 0);

	set<int>driver_set;


	vector<int> left_alternative;
	vector<int>right_alternative;

	for (int ix = 1; ix <= sample_times; ++ix) {
		init_target_node_map_to_original();
		control_link.clear();
		update_control_link();

		if (ix % 500 == 0) {
			printf("第%d次采样!\n", ix);
		}

		if (left_remain.size() > 0) {
			int pick_1 = rand() % left_remain.size();
			int left_pick_node = left_remain[pick_1];
			//cout << "left picked node:" << left_pick_node << endl;
			left_alternative = find_alternative_node_set(LEFT, left_pick_node, cur_match_left);

			//display_1(left_alternative);
			int new_left_matched_node = left_alternative[rand() % left_alternative.size()];
			left_remain[pick_1] = new_left_matched_node;
			unmatch_removed_node_and_max_matching_new_node(LEFT, left_pick_node, new_left_matched_node, cur_match_left, cur_match_right);
		}
		if (right_remain.size() > 0) {
			int pick_2 = rand() % right_remain.size();
			int right_pick_node = right_remain[pick_2];


			right_alternative = find_alternative_node_set(RIGHT, right_pick_node, cur_match_right);
			int new_right_matched_node = right_alternative[rand() % right_alternative.size()];
			right_remain[pick_2] = new_right_matched_node;

			//cout << "right picked node:" << right_pick_node << "\t" << new_right_matched_node << endl;
			unmatch_removed_node_and_max_matching_new_node(RIGHT, right_pick_node, new_right_matched_node, cur_match_left, cur_match_right);
		}

		left_match_set = merge_remain_with_always_matched(left_remain, left_always);
		right_match_set = merge_remain_with_always_matched(right_remain, right_always);
		driver_set = get_right_unmatched_nodes(target_set, right_match_set);

		update_target_node_map(target_set, cur_match_right);
		set<int>duplicate_target_set = update_control_link();
		remove_repeated_nodes(left_match_set, duplicate_target_set);
		random_max_matching(left_match_set, driver_set, link_table);
		in_link = link_table;


		for (int x : driver_set) {
			++matching_count[x];
		}
	}
	return matching_count;
}

std::vector<int> Graph::full_contrl_node_sampling(int sample_times)
{
	printf("%d sampling is beginning\n", sample_times);
	set<int>new_target_set;
	for (int i = 0; i < num_node; ++i) {
		new_target_set.insert(i);
	}

	pair<set<int>, set<int>> out_in_match_set = max_matching(new_target_set);
	set<int> right_match_set = out_in_match_set.second;

	vector<int> cur_match_left = match_left;
	vector<int> cur_match_right = match_right;

	set<int>right_always = find_right_always_matched_nodes(cur_match_right, right_match_set, new_target_set);
	vector<int>right_remain = remove_always_matched_nodes_from_matched_set(right_match_set, right_always);

	vector<int>matching_count(num_node, 0);
	vector<int>alternative_set;
	set<int>driver_set;
	for (int ix = 1; ix <= sample_times; ++ix) {
		if (ix % 500 == 0) {
			printf("第%d次采样!\n", ix);
		}
		if (right_remain.size() > 0) {
			int pick = rand() % right_remain.size();
			int pick_node_1 = right_remain[pick];
			alternative_set = find_alternative_node_set(RIGHT, pick_node_1, cur_match_right);
			int new_right_matched_node = alternative_set[rand() % alternative_set.size()];
			right_remain[pick] = new_right_matched_node;
			//cout << pick_node_1 << "\t" << new_right_matched_node << endl;
			unmatch_removed_node_and_max_matching_new_node(RIGHT, pick_node_1, new_right_matched_node, cur_match_left, cur_match_right);
		}

		right_match_set = merge_remain_with_always_matched(right_remain, right_always);
		driver_set = get_right_unmatched_nodes(new_target_set, right_match_set);

		for (int x : driver_set) {
			++matching_count[x];
		}
	}
	return matching_count;
}
