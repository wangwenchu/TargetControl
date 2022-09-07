#pragma once

#include "Util.h"
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <thread>

using namespace std;

enum SIDE { LEFT = 0, RIGHT };


class Graph {
private:
	int num_node;
	int num_edge;
	vector<vector<int>> in_link;
	vector<vector<int>> out_link;
	vector<int>mark;
	vector<int> match_left;
	vector<int> match_right;
	set<int>target_set;
	unordered_map<int, int>target_node_map;
	unordered_map<int, unordered_set<int>> control_link;
public:
	Graph();
	void read_edgelist(const string& gpath);
	void read_target_set(const string& target_set_path);
	void random_select_target_node(double fraction);
public:
	int get_num_of_nodes()const { return num_node; }
	int get_target_set_size() const { return target_set.size(); }
	vector<vector<int>> get_inLink() const { return in_link; }
	vector<vector<int>> get_outLink() const { return out_link; }
	vector<int> get_match_left() const { return match_left; }
	vector<int> get_match_right() const { return match_right; }
	vector<int> target_control_node_sampling(int sample_times);
	vector<int> full_contrl_node_sampling(int sample_times);
private:
	pair<set<int>, set<int>> max_matching(const set<int>& new_target_set);
	pair<set<int>, set<int>> full_control();
	set<int> find_left_always_matched_nodes(const vector<int>& cur_match_list, const set<int>& match_set);
	set<int> find_right_always_matched_nodes(const vector<int>& cur_match_list, const set<int>& match_set, const set<int>& cur_target_set);
	set<int> get_right_unmatched_nodes(const set<int>& cur_target_set, const set<int>& right_match_set);
	vector<int>remove_always_matched_nodes_from_matched_set(const set<int>& matched_set, const set<int>& always_matched);
	vector<int> find_alternative_node_set(int select_side, int target_node, const vector<int>& cur_match_list);
	void alloc_memory(size_t graphsize);
	bool left_dfs(int u);
	bool right_dfs(int u);
	vector<int>find_left_altertive_set(int target_node, const vector<int>& cur_match_list);
	vector<int>find_right_altertive_set(int target_node, const vector<int>& cur_match_list);
	set<int>merge_remain_with_always_matched(const vector<int>& remain_set, const set<int>& alway_match);
	bool isSubset(const set<int>& small_set, const set<int>& big_set);
	void remove_repeated_nodes(set<int>& s1, const set<int>& s2);
	void init_target_node_map_to_original();
	void update_target_node_map(const set<int>& cur_target_set, const vector<int>& cur_match_right);
	set<int> update_control_link();
	vector<vector<int>> randomize_node_label_and_link_table(const vector<vector<int>>& link_table);
	set<int> get_maped_node(const set<int>& s, const vector<int>& node_pool);
	vector<int> get_maped_node(const vector<int>& s, const vector<int>& node_pool);
	void unmatch_removed_node_and_max_matching_new_node(int select_side, int unmatch_node, int new_pick_node, vector<int>& cur_match_left, vector<int>& cur_match_right, const set<int>& always_match_set, const vector<int>& remain_set);
	void unmatch_left_node_and_max_matching(int unmatch_node, int new_pick_node, vector<int>& cur_match_left, vector<int>& cur_match_right, const set<int>& always_match_set, const vector<int>& remain_set);
	void unmatch_right_node_and_max_matching(int unmatch_node, int new_pick_node, vector<int>& cur_match_left, vector<int>& cur_match_right, const set<int>& always_match_set, const vector<int>& remain_set);
	void mark_for_max_matching_a_node(const set<int>& always_match_set, const vector<int>& remain_set);
	void random_max_matching(set<int>& new_target_set, set<int>& driver_set, const vector<vector<int>>& link_table);
};

