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

struct
{
	vector<vector<int>> link_to_left;
	vector<vector<int>> link_to_right;
	vector<int> temp_match_right;       //temp_match_right[1] = 3:     3(out) -- 1(in) is matched 
	vector<int> temp_match_left;        //temp_match_left[3] = 1:     3(out) -- 1(in) is matched 
	vector<int> match_right;         //similar to temp_match_right
	vector<int> match_left;        //similar to temp_match_left
	vector<int> right_mark;        //mark nodes in in-set, right_mark[4] = 1: node 4 is visited
	vector<int> left_mark;
	int number_of_node = 0;    //网络节点数
	set<int> original_target;   //original target set
	set<int> new_target;       //current target set
	int max_node = 0;      // max id of nodes
	vector<int>all_nodes;
	map<int, set<int>>target_link; //control path
	map<int, int>new_target_to_original; // new_target_to_original[u] = v: v is the orginal target node in  control path of node u
} link;



template <class T>
void init(T& v) 
{
	for (int i = 0; i < int(v.size()); ++i)
	{
		v[i] = 0;
	}
}

int getMaxID(const string& filename) {
	ifstream fin(filename);
	int a, b;
	int maxid = 0;
	while (fin >> a >> b) {
		maxid = max(maxid, max(a, b));
		if (a == 0 || b == 0) {
			cerr << " edge error!" << endl;
			exit(1);
		}
	}
	fin.close();
	return maxid;
}

void read_edgelist(string filename, int size)
{
	link.link_to_left.resize(size);
	link.link_to_right.resize(size);

	ifstream fin1(filename);
	int a, b;
	int number_of_edge = 0;
	while (fin1 >> a >> b)
	{
		link.link_to_left[a].push_back(b);
		link.link_to_right[b].push_back(a);
		link.number_of_node = max(link.number_of_node, max(a, b));
		++number_of_edge;
	}
	fin1.close();

	for (int i = 0; i < size; ++i) {
		link.all_nodes.push_back(i);
	}
	cout << "number of nodes:" << link.number_of_node << endl;
	cout << "number of links:" << number_of_edge << endl;
}

void read_target_node(string filename)
{
	ifstream fin2(filename);
	int a;
	while (fin2 >> a)
	{
		link.original_target.insert(a);
		link.new_target.insert(a);
	}
	fin2.close();
	cout << "number of target node:" << link.new_target.size() << endl;
}

void init_target_map_to_self(const set<int>& target) {
	link.new_target_to_original.clear();
	for (auto x : target) {
		link.new_target_to_original[x] = x;
	}
}

set<int> get_all_neighbor_and_reachable_nodes_to_current_target(const set<int> new_target)
{
	set<int> temp;
	for (auto t : new_target)
	{
		temp.insert(link.link_to_right[t].begin(), link.link_to_right[t].end());
	}
	return temp;
}

/*
* start from node u in out-set to find an augumenting path 
*/
bool Left_DFS(int u) 
{

	for (int i = 0; i < int(link.link_to_left[u].size()); ++i)
	{
		int t = link.link_to_left[u][i];
		if (link.right_mark[t] == false)
		{
			link.right_mark[t] = true;
			if (link.temp_match_right[t] == 0 || Left_DFS(link.temp_match_right[t]) == true)
			{
				link.temp_match_right[t] = u;
				link.temp_match_left[u] = t;
				return true;
			}
		}
	}
	return false;
}

/*
* start from node u in in-set to find an augumenting path
*/
bool Right_DFS(int u) 
{
	for (int i = 0; i < int(link.link_to_right[u].size()); ++i)
	{
		int v = link.link_to_right[u][i];
		if (link.left_mark[v] == false)
		{
			link.left_mark[v] = true;
			if (link.temp_match_left[v] == 0 || Right_DFS(link.temp_match_left[v]) == true)
			{
				link.temp_match_left[v] = u;
				link.temp_match_right[u] = v;
				return true;
			}
		}
	}
	return false;
}

set<int> full_control()
{
	set<int> right_unmatch;
	init(link.temp_match_left);
	init(link.temp_match_right);
	for (int i = 1; i <= link.max_node; ++i)
	{
		std::fill(link.right_mark.begin(), link.right_mark.end(), 0);
		Left_DFS(i);
	}
	for (int j = 1; j <= link.max_node; ++j)
	{
		if (link.temp_match_right[j] == 0)
		{
			right_unmatch.insert(j);
		}
	}
	return right_unmatch;
}

int find_cacti(set<int> right_unmatch)
{
	int size = right_unmatch.size();
	vector<vector<int>> cacti_link(size);
	int u;
	int index = 0;
	int count = 0;
	for (auto t : right_unmatch)
	{
		u = t;
		count++;
		cacti_link[index].push_back(u);
		while (link.temp_match_left[u] != 0)
		{
			cacti_link[index].push_back(link.temp_match_left[u]);
			u = link.temp_match_left[u];
		}
		++index;
	}
	set<int>s;
	int num_of_cactic = 0;
	for (auto t : cacti_link)
	{
		for (auto k : t)
		{
			if (link.new_target.count(k) == 1)
			{
				++num_of_cactic;
				break;
			}

		}
	}

	for (auto t : cacti_link[0])
	{
		s.insert(t);
	}

	return num_of_cactic;
}


/*
* return matched sets in out-set and in-set
*/
pair<set<int>, set<int>> max_matching(const set<int> &left_to_target, const set<int> &new_target)
{
	init(link.temp_match_left);
	init(link.temp_match_right);
	init(link.match_left);
	init(link.match_right);

	set<int> left_match_set;  
	set<int> right_match_set; 
	pair<set<int>, set<int>> left_right_match;

	/*
	* mark all nodes to forbiden other nodes to visit them
	*/
	std::fill(link.right_mark.begin(), link.right_mark.end(), 1);   
	std::fill(link.left_mark.begin(), link.left_mark.end(), 1);
	for (auto t : new_target) {
		for (auto r : left_to_target) {
			link.left_mark[r] = false;
		}
		if (Right_DFS(t)) {
			right_match_set.insert(t);
		}
	}
	for (auto t : right_match_set) {
		left_match_set.insert(link.temp_match_right[t]);
	}

	copy(begin(link.temp_match_left), end(link.temp_match_left), begin(link.match_left));
	copy(begin(link.temp_match_right), end(link.temp_match_right), begin(link.match_right));


	left_right_match.first = left_match_set; //左边的在前
	left_right_match.second = right_match_set;
	return left_right_match; //右边匹配的节点
}



/*
finding always matched nodes in in-set
 */

set<int> finding_always_matching_node_on_the_right(const set<int> &right_match)
{
	set<int> right_always;
	for (auto t : right_match)
	{
		for (int i = 1; i < link.max_node; ++i)
		{
			link.temp_match_right[i] = link.match_right[i];
		}
		for (auto k : link.new_target) //其他节点在最大匹配时已经被标记为true了，而为true的不会变为false
		{
			link.right_mark[k] = false;
		}
		link.right_mark[t] = true;
		if (Left_DFS(link.match_right[t]) == false)
		{
			right_always.insert(t);
		}
	}
	return right_always;
}

/*
finding always matched nodes in out-set
  */
set<int> finding_always_matching_node_on_the_left(const set<int> &left_match)
{
	set<int> left_always;
	for (auto t : left_match)
	{
		//copy(begin(link.match_left), end(link.match_left), begin(link.temp_match_left));
		for (int i = 1; i <= link.max_node; ++i)
		{
			link.temp_match_left[i] = link.match_left[i];
		}
		//memset(link.left_mark, 0, sizeof(link.left_mark));
		std::fill(link.left_mark.begin(), link.left_mark.end(), 0);
		link.left_mark[t] = true;
		if (Right_DFS(link.match_left[t]) == false)
		{
			left_always.insert(t);
		}
	}
	return left_always;
}

set<int> finding_unmatched_node(const set<int>& right_match_set) //找到每次迭代右边不匹配的
{
	set<int> unmatched_set;
	set_difference(link.new_target.begin(), link.new_target.end(), right_match_set.begin(), right_match_set.end(), inserter(unmatched_set, unmatched_set.begin()));
	return unmatched_set;
}

vector<int> remove_always_matching_node_from_matching_set(const set<int>& matching_set, const set<int>& always_matching)
{
	vector<int> remain_node;
	for (auto t : matching_set)
	{
		if (always_matching.count(t) == 0)
		{
			remain_node.push_back(t);
		}
	}
	return remain_node;
}

/*
* node u is to be replaced, here we find the alternative node set that can replace u
*/
pair<vector<int>, vector<pair<vector<int>, vector<int>>>> finding_alternative_node_of_right_selected_node(int u) //这个节点必须满足在taget node中
{
	
	vector<int> alternative;
	int remove = u;
	int node = link.match_right[u];

	vector<int> last_match(link.max_node + 1);
	copy(begin(link.match_right), end(link.match_right), begin(link.temp_match_right));
	copy(begin(link.match_right), end(link.match_right), begin(last_match));

	link.temp_match_right[remove] = 0;

	//only nodes in new target can be visited
	std::fill(link.right_mark.begin(), link.right_mark.end(), 1);
	for (auto it : link.new_target)
	{
		link.right_mark[it] = false;
	}
	link.right_mark[remove] = true;  //making node u inaccessible

	vector<pair<vector<int>, vector<int>>> crp_left_right_match; 
	vector<int> crp_match_left(link.max_node + 1, 0);
	vector<int> crp_match_right(link.max_node + 1, 0);

	pair<vector<int>, vector<int>> crp_p;
	int flag;
	while (Left_DFS(node) == true)
	{
		flag = 0;
		for (auto t : link.new_target)
		{
			if (link.temp_match_right[t] == 0)
			{
				continue;
			}
			crp_match_left[link.temp_match_right[t]] = t;
			if (flag == 1)
				continue;
			if (last_match[t] == 0) 
			{
				//cout << "remove：" << t << endl;
				alternative.push_back(t);
				remove = t;                          
				node = link.temp_match_right[remove]; 
				++flag;
			}
		}
		if (flag != 1)
		{
			cout << "error in finding right alternative 1" << endl;
		}
		crp_match_right = link.temp_match_right;
		crp_p.first = crp_match_left;
		crp_p.second = crp_match_right;
		crp_left_right_match.push_back(crp_p); //crp_left_right_match[index]与节点alternative[index]对应，表示当时连边的匹配数组

		copy(begin(link.temp_match_right), end(link.temp_match_right), begin(last_match));
		link.temp_match_right[remove] = 0; //让新的节点不匹配

		
		std::fill(link.right_mark.begin(), link.right_mark.end(), 1);
		for (auto ix : link.new_target) 
		{
			link.right_mark[ix] = false;
		}
		for (auto k : alternative) 
		{
			link.right_mark[k] = true;
		}
		link.right_mark[u] = true;
	}

	if (alternative.size() == 0)
	{
		cout << "error in finding right alternative 2!" << endl;
	}

	pair<vector<int>, vector<pair<vector<int>, vector<int>>>> alternative_match_arr;
	alternative_match_arr.first = alternative;
	alternative_match_arr.second = crp_left_right_match;
	return alternative_match_arr;
}


pair<vector<int>, vector<pair<vector<int>, vector<int>>>> finding_alternative_node_of_left_selected_node(int u, set<int> left_to_new_target)
{
	vector<int> alternative;

	int remove = u;
	int node = link.match_left[u];

	//int last_match[maxm];
	vector<int>last_match(link.max_node + 1);
	copy(begin(link.match_left), end(link.match_left), begin(link.temp_match_left));
	copy(begin(link.match_left), end(link.match_left), begin(last_match));
	link.temp_match_left[remove] = 0;

	//fill(link.left_mark, link.left_mark + maxm, false);
	std::fill(link.left_mark.begin(), link.left_mark.end(), 0);
	link.left_mark[u] = true;

	vector<pair<vector<int>, vector<int>>> crp_left_right_match; //每个alternative对应的左右匹配标记数组
	vector<int> crp_match_left(link.max_node + 1, 0);
	vector<int> crp_match_right(link.max_node + 1, 0);

	pair<vector<int>, vector<int>> crp_p;
	int flag = 0;
	while (Right_DFS(node) == true)
	{
		flag = 0;
		for (auto t : left_to_new_target)
		{
			if (link.temp_match_left[t] == 0)
				continue;
			crp_match_right[link.temp_match_left[t]] = t;
			if (flag == 1)
				continue;
			if (last_match[t] == 0) //原来没有匹配，现在匹配了
			{
				alternative.push_back(t);
				remove = t;                          //
				node = link.temp_match_left[remove]; //新的增广路径起点
				++flag;
			}
		}
		if (flag != 1)
		{
			cout << "error in finding left alternative 1" << endl;
		}
		crp_match_left = link.temp_match_left;
		crp_p.first = crp_match_left;
		crp_p.second = crp_match_right;
		crp_left_right_match.push_back(crp_p);

		copy(begin(link.temp_match_left), end(link.temp_match_left), begin(last_match));
		link.temp_match_left[remove] = 0;
		//fill(link.left_mark, link.left_mark + maxm, false);
		std::fill(link.left_mark.begin(), link.left_mark.end(), 0);
		for (auto k : alternative)
		{
			link.left_mark[k] = true;
		}
		link.left_mark[u] = true;
	}
	pair<vector<int>, vector<pair<vector<int>, vector<int>>>> alternative_match_arr;
	alternative_match_arr.first = alternative;
	alternative_match_arr.second = crp_left_right_match;

	return alternative_match_arr;
}


set<int> merge_remain_set_with_always_mathing_set(vector<int>& remain, set<int>& always)
{
	set<int> matching_set;
	matching_set = always;
	for (auto& t : remain)
	{
		matching_set.insert(t);
	}
	return matching_set; 
}

set<int>find_upper_bound_in_one_match(set<int>right_unmatch) { 
	int size = right_unmatch.size();
	vector<vector<int>> cacti_link(size);
	int u;
	int index = 0;
	for (auto t : right_unmatch)
	{
		u = t;
		cacti_link[index].push_back(u);
		while (link.match_left[u] != 0) //
		{
			cacti_link[index].push_back(link.match_left[u]);
			u = link.match_left[u];
		}
		++index;
	}
	set<int>s;
	for (auto t : cacti_link)
	{
		for (auto k : t)
		{
			if (link.original_target.count(k) == 1)
			{
				s.insert(t[0]);
				break;
			}

		}
	}
	return s; //可能为空
}

vector<set<int>> get_upper_bound_in_full_control_sampling(int sampling_times) {
	set<int>match_set; //右边匹配的集合
	set<int>always_match;
	vector<int>right_remain; //去除始终匹配集合后剩余的节点
	set<int>unmatch; //右边不匹配的节点，用它寻找路径源头
	set<int>src;
	vector<set<int>> res;
	int pick1, pick_node1;
	int pick2, pick_node2;

	pair<vector<int>, vector<pair<vector<int>, vector<int>>>>right_alternative_and_left_right_match_arr;
	vector<pair<vector<int>, vector<int>>> left_right_match_arr;
	vector<int>right_alternative;

	init(link.temp_match_left);
	init(link.temp_match_right);
	for (int i = 1; i <= link.max_node; ++i)
	{
		//memset(link.right_mark, 0, sizeof(link.right_mark));
		std::fill(link.right_mark.begin(), link.right_mark.end(), 0);
		Left_DFS(i);
	}
	for (int j = 1; j <= link.max_node; ++j)
	{
		if (link.temp_match_right[j] != 0)
		{
			match_set.insert(j);
		}
		else {
			unmatch.insert(j);
		}
	}
	link.match_right = link.temp_match_right; //先复制匹配数组
	link.match_left = link.temp_match_left;

	always_match = finding_always_matching_node_on_the_right(match_set);
	right_remain = remove_always_matching_node_from_matching_set(match_set, always_match);
	//cout << "always match:" << always_match.size() << "\t match_set" << match_set.size() <<"in up"<< endl;

	if (right_remain.size() == 0) { //最大匹配右边节点都始终被匹配，那么修改不了
		src = find_upper_bound_in_one_match(unmatch);
		res.push_back(src);
		return res;
	}
	link.new_target.clear();
	for (int i = 1; i <= link.max_node; ++i) {
		link.new_target.insert(i);
	}//记住这里link.new_target变了
	src.clear();
	res.clear();

	for (int i = 0; i < sampling_times; ++i) { //能运行这表明remian元素个数大于0
		pick1 = rand() % right_remain.size();
		pick_node1 = right_remain[pick1];

		right_alternative_and_left_right_match_arr = finding_alternative_node_of_right_selected_node(pick_node1);
		right_alternative = right_alternative_and_left_right_match_arr.first;
		left_right_match_arr = right_alternative_and_left_right_match_arr.second;

		assert(right_alternative.size() > 0);
		//if (right_alternative.size() > 0)
		{
			pick2 = rand() % right_alternative.size();
			pick_node2 = right_alternative[pick2];
			right_remain[pick1] = pick_node2;

			//修改匹配数组
			link.match_left = left_right_match_arr[pick2].first;
			link.match_right = left_right_match_arr[pick2].second;

			//找源头节点up
			assert(pick_node1 != pick_node2);
			/*	if (pick_node2 == pick_node1) {
					cout << "error2 in full control sampling！" << endl;
				}*/
		}
		unmatch.erase(pick_node2);
		unmatch.insert(pick_node1);
		src = find_upper_bound_in_one_match(unmatch);
		res.push_back(src);
	}
	return res;
}

//randomize the labels of nodes
pair<map<int, int>, map<int, int>> Randomize_Label_of_nodes(const vector<vector<int>>& RLinkTable) {
	//srand(unsigned(time(NULL)));
	vector<int>temp = link.all_nodes;
	random_shuffle(temp.begin() + 1, temp.end()); //index->value

	pair<map<int, int>, map<int, int>> p;
	map<int, int>mp1;
	map<int, int>mp2;

	for (int i = 1; i <= link.max_node; ++i) {    //original: index,new :value
		mp1[i] = temp[i]; //index->val
		mp2[temp[i]] = i;  //vaule->index  
		//cout << i << "\t" << temp[i] << endl;
	}
	p.first = mp1;
	p.second = mp2;
	
	link.link_to_right.clear();
	link.link_to_right.resize(link.max_node + 1);
	for (int i = 1; i <= link.max_node; ++i) {
		for (auto t : RLinkTable[i]) {
			link.link_to_right[mp1[i]].push_back(mp1[t]);  //t1->i,t2->i,t3->i
		}
		random_shuffle(link.link_to_right[mp1[i]].begin(), link.link_to_right[mp1[i]].end());
	}
	return p;
}


set<int>get_newLabeled_target_nodes(set<int>target_set, map<int, int>mp) {
	set<int> temp;
	for (auto t : target_set) {
		int newnode = mp[t];   //t: orignal,value:new label
		temp.insert(newnode);
	}
	return temp;
}

set<int>get_original_nodeLable(const set<int>& s, map<int, int>& mp) {
	set<int>temp;
	for (auto t : s) {
		int node = mp[t];
		temp.insert(node);
	}
	return temp;
}


void duplicate_nodes_in_driver_set(set<int>& s1, set<int>s2) {
	for (auto x : s2) {
		if (s1.count(x)) {
			s1.erase(x);
		}
	}
}

//判断s1是不是s2的子集
bool IsSubset(const set<int>&s1, const set<int> &s2) {
	assert(!s1.empty() && !s2.empty());
	for (auto x : s1) {
		if (s2.count(x) == 0) {
			return false;
		}
	}
	return  true;
}



void  remove_repeated_nodes_in_new_target(set<int>& s1, const set<int>& s2) {
	for (auto x : s2) {
		assert(s1.count(x) == 1);
		s1.erase(x);
	}
}

vector<int> get_original_match_left_arr(vector<int>match_left, map<int, int>mp) {
	vector<int>res(match_left.size());
	for (int i = 1; i <= link.max_node; ++i) {
		if (match_left[i] > 0) {
			res[mp[i]] = mp[match_left[i]];
		}
	}
	return res;
}



void update_map(const set<int>& left_match, const vector<int>& match_left) {
	map<int, int>mp;
	for (auto t : left_match) {

		int u = match_left[t];   //    v-->t-->u--u
		//     t-->u
		mp[t] = link.new_target_to_original[u];  //t ->u
	}
	link.new_target_to_original = mp;
}

set<int>get_duplicated_nodes_in_control_link(const set<int>& left_match) {
	set<int>has_before;
	for (auto x : left_match) {
		assert(link.new_target_to_original.count(x));
		int u = link.new_target_to_original[x];  //original target
		if (link.target_link[u].count(x)) {
			has_before.insert(x);
			link.target_link.erase(u); //
		}
	}
	return has_before;
}

void update_control_link() {
	for (auto t : link.new_target_to_original) {
		int u = t.second;
		link.target_link[u].insert(t.first);
	}
}



//randomize labels stratrgy
set<int>random_match_inLaterInteration(const vector<vector<int>>& LinkTable, set<int>left_match_set, set<int> first_unamtch) {
	
	set<int>right_match;
	set<int>left_match;
	set<int>right_unmatch;
	pair<set<int>, set<int>>left_right_match;

	set<int>driver_set = first_unamtch;

	set<int>left_to_target;
	link.new_target = left_match_set;

	set<set<int>> mem;
	set<int>last_target;
	pair<map<int, int>, map<int, int>>p;

	while (!link.new_target.empty()) {  

		last_target = link.new_target;
		p = Randomize_Label_of_nodes(LinkTable);
		link.new_target = get_newLabeled_target_nodes(link.new_target, p.first);

		left_to_target = get_all_neighbor_and_reachable_nodes_to_current_target(link.new_target); //正是通过link.link_to_right计算,得到的节点是重新编号的
		left_right_match = max_matching(left_to_target, link.new_target);
		left_match = left_right_match.first;

		right_match = left_right_match.second;
		right_unmatch = finding_unmatched_node(right_match);  //new_target - right_match

		right_unmatch = get_original_nodeLable(right_unmatch, p.second);

		driver_set.insert(right_unmatch.begin(), right_unmatch.end());

		
		left_match = get_original_nodeLable(left_match, p.second);

		if (left_match.size() == last_target.size()) {
			mem.insert(last_target);
		}
		else {
			mem.clear();
		}
	
		if (mem.count(left_match)) {
			link.new_target.clear();
			break;
		}

		link.match_left = get_original_match_left_arr(link.match_left, p.second);

		update_map(left_match, link.match_left); //将左边匹配的节点和最初target节点映射
		set<int>has_before = get_duplicated_nodes_in_control_link(left_match);
		update_control_link();  //control_link:从target到当前迭代的链上的节点，map<int,set>
		
		remove_repeated_nodes_in_new_target(left_match, has_before); //left_match is changed
	
		link.new_target = left_match;

	}
	return driver_set;
}


void get_control_link_nodes(set<int>last_target, vector<int> match_right) {
	for (int x : last_target) {
		int u = match_right[x];
		if (u > 0) {
			link.target_link[x].insert(u);
		}
	}
}

void Add_original_target_nodes() {
	for (auto x : link.original_target) {
		link.target_link[x].insert(x);
	}
}







vector<set<int>> sampling_of_target_control(int sampling_times)
{
	cout << sampling_times << " sampling is beginning" << endl;
	link.new_target = link.original_target;
	set<int>left_to_new_target = get_all_neighbor_and_reachable_nodes_to_current_target(link.new_target);

	pair<set<int>, set<int>>left_right_match = max_matching(left_to_new_target, link.new_target); //第一轮迭代最大匹配
	set<int>left_match_set = left_right_match.first;   //可能为空集
	set<int>right_match_set = left_right_match.second;
	set<int>left_always = finding_always_matching_node_on_the_left(left_match_set);  //可能为空集
	set<int>right_always = finding_always_matching_node_on_the_right(right_match_set);

	vector<int>left_remain = remove_always_matching_node_from_matching_set(left_match_set, left_always);
	vector<int>right_remain = remove_always_matching_node_from_matching_set(right_match_set, right_always);


	//先记录第一轮最大匹配的相关信息

	set<int>firstLMS = left_match_set;
	set<int>firstRMS = right_match_set;
	set<int>first_L_Always = left_always;
	set<int>first_R_Always = right_always;
	vector<int>first_L_remain = left_remain;
	vector<int>first_R_remain = right_remain;
	vector<int>firstML = link.match_left;
	vector<int>firstMR = link.match_right;
	set<int>first_L2Tar = left_to_new_target;


	pair<vector<int>, vector<pair<vector<int>, vector<int>>>>R_alternative_and_LRM_arr;
	pair<vector<int>, vector<pair<vector<int>, vector<int>>>> L_alternative_and_LRM_arr;
	vector<pair<vector<int>, vector<int>>> LRM_arr;

	//后续将打乱编号,邻接表link.link_to_right会被修改,先记录下来
	const vector<vector<int>>LinkTable = link.link_to_right;


	vector<set<int>> list_driver_set;
	set<int>driver_set;
	for (int iter = 0; iter < sampling_times; ++iter) {
		init_target_map_to_self(link.original_target);
		link.target_link.clear(); //清空控制链,并将初始控制的节点加入

		Add_original_target_nodes();

		link.link_to_right = LinkTable;
		if ((iter + 1) % 500 == 0)
		{
			printf("sampled %d times\n", iter + 1);
		}

		link.match_left = firstML;
		link.match_right = firstMR;
		right_remain = first_R_remain;
		left_remain = first_L_remain;
		link.new_target = link.original_target;
		if (right_remain.size()) {
			int pick1 = rand() % right_remain.size();
			int pick_node1 = right_remain[pick1];
			R_alternative_and_LRM_arr = finding_alternative_node_of_right_selected_node(pick_node1);
			vector<int>right_alternative = R_alternative_and_LRM_arr.first;
			LRM_arr = R_alternative_and_LRM_arr.second;
			assert(right_alternative.size() > 0);

			int pick2 = rand() % right_alternative.size();
			int pick_node2 = right_alternative[pick2];

			right_remain[pick1] = pick_node2;

			//修改匹配数组
			link.match_left = LRM_arr[pick2].first;
			link.match_right = LRM_arr[pick2].second;

		}
		if (left_remain.size()) {
			int pick3 = rand() % left_remain.size();
			int pick_node3 = left_remain[pick3];

			L_alternative_and_LRM_arr = finding_alternative_node_of_left_selected_node(pick_node3, first_L2Tar);
			vector<int>left_alternative = L_alternative_and_LRM_arr.first;

			LRM_arr = L_alternative_and_LRM_arr.second;

			assert(left_alternative.size() > 0);
			int pick4 = rand() % left_alternative.size();
			int pick_node4 = left_alternative[pick4];
			//cout << "pick4:" << pick_node4 << endl;
			left_remain[pick3] = pick_node4;
		

			link.match_left = LRM_arr[pick4].first;
			link.match_right = LRM_arr[pick4].second;
		}
		
		//record previous matching
		firstML = link.match_left;
		firstMR = link.match_right;
		first_R_remain = right_remain;
		first_L_remain = left_remain;


		right_match_set = merge_remain_set_with_always_mathing_set(right_remain, first_R_Always);
		//print(right_match_set);
		set<int>first_unmatch = finding_unmatched_node(right_match_set); //new_target - right_match_set

		
		left_match_set = merge_remain_set_with_always_mathing_set(left_remain, first_L_Always);
		//获得控制链上的节点,每次在迭代最大匹配前将target set加入控制链,因为第一轮已经将匹配节点修改了,所以这里映射的节点也需要修改
		//print(left_match_set);
		update_map(left_match_set, link.match_left);
		update_control_link();


		driver_set.clear();
		driver_set = random_match_inLaterInteration(LinkTable, left_match_set, first_unmatch); //后续迭代

		//求first_unmatch和后续不匹配节点的并集
		if (driver_set.size()) {
			list_driver_set.push_back(driver_set);
		}
		//cout << "driver_set:" << driver_set.size() << endl;
	}
	link.link_to_right = LinkTable; //将邻接表修改回来
	return list_driver_set;
}


void Compute_frequencyOfMDNS(string filename, const vector<set<int>> list_driver_set) {
	ofstream fout(filename);
	vector<int> count_times(link.max_node + 1, 0);
	for (auto t : list_driver_set) {
		for (auto k : t) {
			++count_times[k];
		}
	}
	double pc = 0.0;
	int N = list_driver_set.size();
	int num[3] = { 0,0,0 };
	double EPS = 1.0e-6;
	for (int i = 1; i <= link.max_node; ++i) {
		pc = double(count_times[i]) / N;
		assert(pc >= 0 && pc <= 1);
		if (pc >= -EPS && pc <= EPS) {
			++num[0];
		}
		else if ((pc - 1) >= -EPS && (pc -1)<= EPS) {
			++num[2];
		}
		else {
			++num[1];
		}
		fout << i << "\t" << link.link_to_left[i].size() << "\t" << link.link_to_right[i].size() << "\t" << pc << endl;
	}
	fout << "-------------------------------------------------------------------------------" << endl;
	fout << "redundant(0),intermittent(0-1),critical(1):" << endl;
	fout << double(num[0]) / link.max_node << "\t" << double(num[1]) / link.max_node << "\t" << double(num[2]) / link.max_node << endl;
}


void full_control_for_three_types_nodes() {
	link.new_target.clear();
	for (int i = 1; i <= link.max_node; ++i) {
		link.new_target.insert(i);
	}
	set<int>left_to_new_target = get_all_neighbor_and_reachable_nodes_to_current_target(link.new_target);
	pair<set<int>, set<int>>left_right_match = max_matching(left_to_new_target, link.new_target); 
	set<int>left_always = finding_always_matching_node_on_the_left(left_right_match.first);
	set<int>right_always = finding_always_matching_node_on_the_right(left_right_match.second);
	//cout << "left right always matched:" << left_always.size() << "\t" << right_always.size() << endl;
	/*for (int x : right_always) {
		cout << x << "\t";
	}
	cout << endl;*/
	set<int>critical;
	for (int i = 1; i <= link.max_node; ++i) {
		if (0 == link.link_to_right[i].size()) {
			critical.insert(i);
		}
	}
	cout << "full control driver set:" << link.new_target.size() - left_right_match.second.size() << endl;
	cout << "0(redudant) \t 0-1(intermimtent) \t 1(critical)" << endl;
	cout << right_always.size() * 1.0 / link.max_node << "\t" << (link.max_node - int(critical.size()) - right_always.size()) * 1.0 / link.max_node << "\t" << critical.size() * 1.0 / link.max_node << endl;
}

set<int> randomly_generate_target_set(int max_node, float f) {
	vector<int>vec(max_node);
	for (int i = 1; i <= max_node; ++i) {
		vec[i - 1] = i;
	}
	random_shuffle(vec.begin(), vec.end());
	int size = f * max_node;
	set<int> res(vec.begin(), vec.begin() + size);
	return res;
}

int main()
{
	/*
	* isTargetControl == true: control a subset,otherwise, control all nodes 
	*/
	bool isTargetControl = true;  
	//isTargetControl = false;

	/*
	* random_pick_target == true, the target set will be selected randomly
	* otherwise, read it from file
	*/

	bool random_pick_target = true;
	

	srand(unsigned(time(NULL)));

	string netfile = "net.txt";  //input file for edgelist
	string outfile = "ans.txt"; //out file


	int graphSize = getMaxID(netfile) + 1;
	read_edgelist(netfile, graphSize);
	link.temp_match_right.resize(graphSize);
	link.temp_match_left.resize(graphSize);
	link.match_left.resize(graphSize);
	link.match_right.resize(graphSize);
	link.left_mark.resize(graphSize);
	link.right_mark.resize(graphSize);
	link.max_node = graphSize - 1;

	if (isTargetControl)
	{
		if (random_pick_target)
		{
			float frac = 0.2;
			link.new_target = randomly_generate_target_set(link.max_node, frac);
			link.original_target = link.new_target;
			cout << "size of target node:" << link.new_target.size() << endl;
		}
		else {
			/*
			* if the target nodes is give, read from this file
			*/
			string target_file = "target.txt";
			read_target_node(target_file);
		}
		
		
		int sampling_times = 1000 * link.original_target.size();

		vector<set<int>> all_driver_set = sampling_of_target_control(sampling_times);
		if (all_driver_set.size()) {
			Compute_frequencyOfMDNS(outfile, all_driver_set);
		}
		else {
			cout << "There is no driver node set!" << endl;
		}
	}
	else
	{
		full_control_for_three_types_nodes();
	}
	std::system("pause");

	return 0;
}


