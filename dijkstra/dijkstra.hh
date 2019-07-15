#pragma once

#include <vector>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>

namespace Elf {

typedef unsigned int uint32;

using std::vector;
using std::priority_queue;
using std::cout;

struct Edge {
	int src;
	int dst;
	int w; // weight
};

inline
bool operator>(const Edge &a, const Edge &b){
	return a.w < b.w;
}

void print_graph(vector<vector<Edge>> adj){
	for(uint32 i = 0; i < adj.size(); i++){
		cout << i << ": ";
		for(auto &edge: adj[i]){
			cout << edge.dst << "(" << edge.w << ") ";
		}
		cout << "\n";
	}
}

void dijkstra(int nNodes = 100){
	vector<vector<Edge>> adjacency(nNodes);

	std::srand(std::time(NULL));

	// Generate edges
	for(int i = 0; i < 2*nNodes; i++){
		Edge e;
		e.src = std::rand() % nNodes;
		e.dst = std::rand() % nNodes;
		e.w   = std::rand();

		// We prohibit self-edges
		if(e.src == e.dst){
			i--;
			continue;
		} else {
			adjacency[e.src].push_back(e);
		}
	}

	// print_graph(adjacency);

	// Generate begin and end nodes
	int beg = std::rand() % nNodes;
	int end = beg;
	while(end == beg) end = std::rand() % nNodes;

	// cout << beg << " -> " << end << "\n";

	// Main data used in the djikstra alg
	vector<bool> marks(nNodes, false);
	priority_queue<Edge, vector<Edge>, std::greater<Edge>> pq;

	// Add the source node
	marks[beg] = true;
	for(Edge &e: adjacency[beg])
		pq.push(e);

	while(marks[end] == false && pq.empty() == false){
		Edge e = pq.top(); pq.pop();
		int visit = e.dst;

		if(marks[visit] == false){
			marks[visit] = true;
			for(Edge &e: adjacency[visit]){
				if(marks[e.dst] == false)
					pq.push(e);
			}
		}
	}

	if(marks[end] == true){
		cout << "Success.\n";
	} else {
		cout << "Fail.\n";
	}
}


} /* namespace Elf */
