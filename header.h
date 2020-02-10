#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include<set>
#include<queue>
#include<string>
#include<math.h>
#include<time.h>
#include<map>
#define epsilon 1e-10
#define VERTEX_MAXNUM  400000
#define EDGE_MAXNUM 500000
#define FACE_MAXNUM 200000
using namespace std;

class Matrix
{
public:
	double data[4][4];
	Matrix()
	{
		for (int i = 0;i < 4;i++)
			for (int j = 0;j < 4;j++)
				data[i][j] = 0;
	}
	Matrix(double a00, double a01, double a02, double a03, double a10, double a11, double a12, double a13,
		double a20, double a21, double a22, double a23, double a30, double a31, double a32, double a33);
};

class Vector_4
{
public:
	double data[4];
	Vector_4() {}
	Vector_4(double a1, double a2, double a3, double a4);
};

Vector_4* solve(Matrix* m, Vector_4* b);
double det(Matrix* m);
Matrix* getQ(int v);
Matrix* add(Matrix* m1, Matrix* m2);

struct Face
{
	int v[3];
	bool valid;
	Face(int a, int b, int c) {
		v[0] = a; v[1] = b; v[2] = c;
	}
	void update(int a, int b, int c) {
		v[0] = a; v[1] = b; v[2] = c;
	}
	bool operator < (const Face& f) const;
};

class Vertex
{
public:
	int id;
	bool valid;
	double coordinate[3];
	set<int> adjacent_v;
	set<int> planes;
	double Q[4][4];
	void calQ();
	void calQ(int a1, int a2);
	void update(int i, double c1, double c2, double c3);
	Vertex(): valid(true) {}
	Vertex(int i, double c1, double c2, double c3);
};

class Edge
{
public:
	int id;
	int v1, v2;
	double cost;
	double vbar[3];
	void update();
	Edge(int a, int b):v1(a), v2(b) {}
	Edge() {}
};

class Heap
{
public:
	struct cmp {
		bool operator() (Edge X, Edge Y) {
			return X.cost > Y.cost;
		}
	};
	priority_queue<Edge, vector<Edge>, cmp> queue;
	map<pair<int, int>, int> ver2edge;
	bool isDeleted[EDGE_MAXNUM];
	int size;
	void push(Edge& e);
	void remove(Edge e);
	void remove(int a, int b);
	Edge pop();
	Heap() :size(0) {}
};
bool isEdge(int v1, int v2);
int getPlaneNum(int v1, int v2);
void readinput(string filename="");
void simplify();
void heapify();
void output(string filename="");