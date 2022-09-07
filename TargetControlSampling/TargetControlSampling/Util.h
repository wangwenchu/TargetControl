#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>


static void compute_frequence_and_write(const std::string& filename, const std::vector<int>& matching_count, int sampling_times) {
	std::ofstream fout(filename);
	if (!fout.good()) {
		printf("output file error!\n");
		std::this_thread::sleep_for(std::chrono::milliseconds(2000)); //µ•Œª «∫¡√Î
		exit(1);
	}
	int n = matching_count.size();
	for (int i = 0; i < n; ++i) {
		fout << i << "\t" << matching_count[i] * 1.0 / sampling_times << std::endl;
	}
	fout.close();
}

template<typename T>
static void display_1(T v) {
	printf("*****************\n");
	for (auto x : v) {
		printf("%d\t", x);
	}
	printf("\n");
	printf("-----------------\n");
}

template<typename T>
static void display_2(T v) {
	for (auto x : v) {
		for (auto t : x) {
			printf("%d\t", x);
		}
		if (x.size() > 0) {
			printf("\n");
		}
	}
}