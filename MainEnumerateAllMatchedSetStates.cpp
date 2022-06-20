#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<set>
#include<string>
#include<algorithm>

using namespace std;

struct Graph {
	vector<vector<int>> link_to_left;
	vector<vector<int>> link_to_right;

	vector<int>mark_left;
	vector<int>match_left;

	vector<int>mark_right;
	vector<int>match_right;
	int n;
	Graph() :n(0) {}
}graph;


struct Target
{
	set<int>s0;
	vector<int> s1; //remove 0 in degree in target
};


int getMaxID(string& filename) {  //label from 0 
	ifstream fin(filename);
	if (!fin.good()) {
		cerr << "Input file doesn't exist or is broken!" << endl;
		exit(0);
	}
	int a, b;
	int maxID = -1;   
	while (fin >> a >> b) {
		maxID = max(maxID, max(a, b));
	}
	fin.close();
	return maxID;
}


void read_edgelist(string& filename) {
	ifstream fin(filename);
	if (!fin.good()) {
		cerr << "Input file doesn't exist or is broken!" << endl;
		exit(0);
	}
	int a, b;
	int m = 0;
	while (fin >> a >> b) {
		graph.link_to_left[a].push_back(b);
		graph.link_to_right[b].push_back(a);
		++m;
	}
	fin.close();
	cout << "number of nodes:" << graph.n << endl;
	cout << "number of edges:" << m << endl;
}


set<int> read_target_node(string& filename) {
	ifstream fin(filename);
	if (!fin.good()) {
		cerr << "Input file doesn't exist or is broken!" << endl;
		exit(0);
	}
	int a;
	set<int>S;
	while (fin >> a) {
		S.insert(a);
	}
	fin.close();
	cout << "number of target nodes:" << S.size() << endl;
	return S;
}

bool Right_DFS(int u) {
	for (int x : graph.link_to_right[u]) {
		if (graph.mark_left[x] == 0) {
			graph.mark_left[x] = 1;
			if (graph.match_left[x] == -1 || Right_DFS(graph.match_left[x])) {
				graph.match_left[x] = u;
				return true;
			}
		}
	}
	return false;
}

bool Left_DFS(int u) {
	for (int x : graph.link_to_left[u]) {
		if (graph.mark_right[x] == 0) {
			graph.mark_right[x] = 1;
			if (graph.match_right[x] == -1 || Left_DFS(graph.match_right[x])) {
				graph.match_right[x] = u;
				return true;
			}
		}
	}
	return false;
}

pair<set<int>,set<int>> maxMatching(const set<int>& target) {
	for (int i = 0; i < graph.match_left.size(); ++i) {
		graph.match_left[i] = -1;
		graph.match_right[i] = -1;
	}
	for (int x : target) {
		for (int i = 0; i < graph.n; ++i) {
			graph.mark_left[i] = 0;
		}
		Right_DFS(x);
	}
	set<int>left_match;
	set<int>right_match;
	for (int i = 0; i < graph.n; ++i) {
		if (graph.match_left[i] != -1) {
			int u = graph.match_left[i];  //i -- u
			graph.match_right[u] = i;

			left_match.insert(i);
			right_match.insert(u);

		}
	}
	auto res = make_pair(left_match, right_match);
	return res;
}


set<int> get_right_always_matching_set(const set<int>&right_match,const vector<int>&match_right,const vector<int>& target) {
	set<int>right_always;
	vector<int>original_match_right = match_right;
	for (int i = 0; i < graph.n; ++i) {
		graph.mark_right[i] = 1;
	}
	for (int x : right_match) {
		graph.match_right = original_match_right;  //keep current matching, and make x inreachable
		graph.match_right[x] = -1;

	
		for (int y :target) {
			graph.mark_right[y] = 0;
		}

		graph.mark_right[x] = 1;
		//cout << "match_right[x]:" << match_right[x]  <<" " << x<< endl;
		if (Left_DFS(original_match_right[x])==false) {
			right_always.insert(x);
		}
	}
	return right_always;
}

set<int>get_left_always_matching_set(const set<int>& left_match, const vector<int> &match_left) {
	set<int>left_always;
	vector<int>original_left_match = match_left;
	for (int x : left_match) {
		graph.match_left = original_left_match;
		graph.match_left[x] = -1;
		int u = match_left[x];

		for (int i = 0; i < graph.n; ++i) {
			graph.mark_left[i] = 0;
		}
		graph.mark_left[x] = 1;
		if (Right_DFS(original_left_match[x]) == false) {
			left_always.insert(x);
		}
	}
	return left_always;
}

vector<int> remove_0_in_degree_in_target(const set<int> &target){
	vector<int>v;
	for (int x : target) {
		if (graph.link_to_right[x].size() != 0) {
			v.push_back(x);
		}
	}
	cout << "not 0 in degree target nodes:" << v.size() << endl;
	return v;
}


vector<vector<int>> construct_set_for_s1_target(const vector<int>&target,const set<int> &left_always) {
	vector<vector<int>>S;

	for (int x: target) {
		vector<int>temp;
		for (int t : graph.link_to_right[x]) {
			if (left_always.count(t)==0) {
				temp.push_back(t);
			}
		}
		if (!temp.empty()) {
			S.push_back(temp);
		}
	}
	return S;
}

vector < vector<int>> construct_set_for_s1_target_without_remove_always(const vector<int>& target) {
	vector<vector<int>>S;
	for (int x : target) {
		S.push_back(graph.link_to_right[x]);
	}
	return S;
}

vector<vector<int>>construct_set_for_left_side(const set<int>& target,const set<int>&right_always) {
	vector<vector<int>>S;
	set<int>remain;
	std::set_difference(target.begin(), target.end(), right_always.begin(), right_always.end(), inserter(remain,remain.begin()));
	for (int i = 0; i < graph.n; ++i) {
		vector<int>temp;
		for (int x : graph.link_to_left[i]) {
			if (remain.count(x)) {
				temp.push_back(x);
			}
		}
		if (!temp.empty()) {
			S.push_back(temp);
		}
	}
	return S;
}


class EnumerateMMS {
private:
	vector<vector<int>>all_set;
public:
	EnumerateMMS(const vector<vector<int>>& v1) :all_set(v1) {}
	void backtraverse(set<set<int>>& MMS, vector<int>& temp,int k, int idx);
};

void EnumerateMMS::backtraverse(set<set<int>>&MMS, vector<int>& temp, int k, int idx) {
	if (temp.size() >= k) {
		set<int>tmp(temp.begin(), temp.end());
		if (tmp.size() == k) {
			MMS.insert(tmp);
			return; //找到了可以直接退出,后面没必要继续找了
		}
	}
	if (idx == all_set.size()) {
		return;
	}
	for (int x : all_set[idx]) {
		temp.push_back(x);
		backtraverse(MMS,temp, k, idx + 1);
		temp.pop_back();
	}
}

class EnumrateStates {
public:
	using SPSS = set < pair<set<int>, set<int>> >;
	vector<int>target;
	vector<vector<int>>setTotarget;
	SPSS all_states;

	void backTraverse(vector<int>& temp, int k,int idx);
	void help(vector<int>& select, set<int>& MMS);
	void DFS(const set<int>&MMS,const vector<set<int>>& S, set<int>temp,int k, int idx);
};

void EnumrateStates::backTraverse(vector<int>& temp, const int k,int idx) {
	
	if (idx == setTotarget.size()) {
		if (temp.size() >= k) {
			set<int>newtemp(temp.begin(), temp.end());  //3,7,4,4->3,4,7
			if (newtemp.size() == k) {
				help(temp, newtemp);
			}
		}
		return;
	}
	for (int t : setTotarget[idx]) {
		temp.push_back(t);
		backTraverse(temp, k, idx + 1);
		temp.pop_back();
	}
}

//select is the nodes selected from out-set

void EnumrateStates::help(vector<int>& select, set<int>& MMS) {
	vector<set<int>>S;
	for (int x : MMS) {
		set<int>temp;
		for (int i = 0; i < select.size(); ++i) {
			if (x == select[i]) {  
				temp.insert(target[i]);
			}
		}
		S.push_back(temp);
	}
	set<int>tmp;
	DFS(MMS, S, tmp, MMS.size(), 0);
}
               

void EnumrateStates::DFS(const set<int>& MMS, const vector<set<int>>& S, set<int>temp, int k, int idx){
	if (temp.size() == k) {
		auto res = make_pair(temp, MMS);
		all_states.insert(res);
	}
	if (idx == S.size()) {
		return;
	}

	for (int t : S[idx]) {
		temp.insert(t);
		DFS(MMS, S, temp, k, idx + 1);
		temp.erase(t);
	}
}


int main()
{
	string filename1 = "net.txt";
	string filename2 = "target.txt";
	
	int graphSize = getMaxID(filename1) + 1;

	graph.n = graphSize;
	graph.link_to_left.resize(graphSize);
	graph.link_to_right.resize(graphSize);
	graph.mark_left.resize(graphSize);
	graph.match_left.resize(graphSize);
	graph.mark_right.resize(graphSize);
	graph.match_right.resize(graphSize);
	read_edgelist(filename1);

	Target tar;
	tar.s0 = read_target_node(filename2); 
	tar.s1 = remove_0_in_degree_in_target(tar.s0);
	
	auto match_set = maxMatching(tar.s0);
	
	set<int>left_match = match_set.first;
	set<int>right_match = match_set.second;
	int matchingSize = right_match.size();

	cout << "maximum matching size:" << matchingSize << endl;
	
	set<int>right_always = get_right_always_matching_set(right_match,graph.match_right,tar.s1);
	cout << "right always matched set:" << right_always.size() << endl;

	set<int>left_always = get_left_always_matching_set(left_match, graph.match_left);
	cout << "left always matched set:" << left_always.size() << endl;

	auto S0 = construct_set_for_left_side(tar.s0,right_always);
	auto S1 = construct_set_for_s1_target(tar.s1,left_always);

	EnumerateMMS obj0(S0);
	set<set<int>>MMS_in;
	vector<int>temp;
	obj0.backtraverse(MMS_in, temp, matchingSize - right_always.size(), 0);


	EnumerateMMS obj1(S1);
	set<set<int>>MMS_out;
	temp.clear();
	obj1.backtraverse(MMS_out, temp, matchingSize - left_always.size(), 0);

	cout << "number of all left matched set:" << MMS_out.size() << endl;
	cout << "number of all right matched set:" << MMS_in.size() << endl;

	EnumrateStates obj;
	obj.target = tar.s1;
	obj.setTotarget = construct_set_for_s1_target_without_remove_always(tar.s1);
	vector<int>tmp;
	obj.backTraverse(tmp, matchingSize, 0);
	cout << "number of all states:" << obj.all_states.size() << endl;
	system("pause");
	return 0;
}