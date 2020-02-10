#include "header.h"

Vertex vertices[VERTEX_MAXNUM];
int vertex_num = 0;
int face_num = 0;
int goal_num = 0;
int start_check = 0;
double rate = 0.03;
vector<Face> faces;
Heap heap;

Vector_4::Vector_4(double a1, double a2, double a3, double a4)
{
	data[0] = a1;	data[1] = a2;	data[2] = a3;	data[3] = a4;
}

Vector_4 * solve(Matrix* m, Vector_4* b)
{
	double a[4][4];
	for (int i = 0;i < 4;i++)
		for (int j = 0;j < 4;j++)
			a[i][j] = m->data[i][j];
	double res[4];
	
	for (int i = 0;i <  3;i++)
	{
		//if a[i][i] == 0
		if (abs(a[i][i]) < epsilon)
		{
			for (int k = i + 1;k < 4;k++)
			{
				if (abs(a[k][i]) > epsilon)
				{
					for (int l = i;l < 4;l++)
					{
						double temp = a[i][l];
						a[i][l] = a[k][l];
						a[k][l] = temp;
					}
					double temp = b->data[k];
					b->data[k] = b->data[i];
					b->data[i] = temp;
					break;
				}
			}
		}
		for (int j = i + 1;j < 4;j++)
		{
			double e = a[j][i] / a[i][i];
			for (int k = i;k < 4;k++)
			{
				a[j][k] -= e * a[i][k];
			}
			b->data[j] -= e * b->data[i];
		}
	}
	for (int i = 3;i > 0;i--)
	{
		for (int j = i - 1;j >= 0;j--)
		{
			double e = a[j][i] / a[i][i];
			a[j][i] = 0;
			b->data[j] -= e * b->data[i];
		}
	}
	for (int i = 0;i < 4;i++)
	{
		res[i] = b->data[i] / a[i][i];
	}

	Vector_4* dest = new Vector_4(res[0], res[1], res[2], res[3]);
	return dest;
}

double det(Matrix * m)
{
	double a[5][5];
	for (int i = 0;i < 4;i++)
		for (int j = 0;j < 4;j++)
			a[i+1][j+1] = m->data[i][j];
	double a1 = a[2][2] * (a[3][3] * a[4][4] - a[3][4] * a[4][3]) - a[3][2] * (a[2][3] * a[4][4] - a[2][4] * a[4][3]) + a[4][2] * (a[2][3] * a[3][4] - a[2][4] * a[3][3]);
	double a2 = a[1][2] * (a[3][3] * a[4][4] - a[3][4] * a[4][3]) - a[3][2] * (a[1][3] * a[4][4] - a[1][4] * a[4][3]) + a[4][2] * (a[1][3] * a[3][4] - a[1][4] * a[3][3]);
	double a3 = a[1][2] * (a[2][3] * a[4][4] - a[2][4] * a[4][3]) - a[2][2] * (a[1][3] * a[4][4] - a[1][4] * a[4][3]) + a[4][2] * (a[1][3] * a[2][4] - a[1][4] * a[2][3]);
	double a4 = a[1][2] * (a[2][3] * a[3][4] - a[2][4] * a[3][3]) - a[2][2] * (a[1][3] * a[3][4] - a[1][4] * a[3][3]) + a[3][2] * (a[1][3] * a[2][4] - a[1][4] * a[2][3]);
	return a1 * a[1][1] - a2 * a[2][1] + a3 * a[3][1] - a4 * a[4][1];
}

Matrix* getQ(int v)
{

	Vertex* v1 = vertices + v;
	Matrix* Q = new Matrix;
	for (int i = 0;i < 4;i++)
		for (int j = 0;j < 4;j++)
			Q->data[i][j] = v1->Q[i][j];
	return Q;
}

Matrix * add(Matrix * m1, Matrix * m2)
{
	Matrix* m = new Matrix;
	for (int i = 0;i < 4;i++)
		for (int j = 0;j < 4;j++)
			m->data[i][j] = m1->data[i][j] + m2->data[i][j];
	return m;
}

bool isEdge(int v1, int v2)
{
	Vertex* p = vertices + v1;
	set<int>::iterator it;
	return p->adjacent_v.count(v2) > 0;
	
}

int getPlaneNum(int v1, int v2)
{
	int count = 0;
	set<int>::iterator it;
	Vertex* p1 = vertices + v1;
	for (it = p1->adjacent_v.begin();it != p1->adjacent_v.end();it++)
		if (isEdge(v2, (*it)))
			count++;

	return count;
}

void readinput(string filename)
{
	FILE *file;
	file = fopen("D:/models/bunny.obj", "r");

	char c = 0;
	int i = 0;
	double x, y, z;
	int a, b, cc;
	do
	{
		c = fgetc(file);
		if (c == '#')
			do
			{
				c = fgetc(file);
			} while (c != '\n');
	} while (c != 'v');

	while (!feof(file))
	{
		if (c != 'v')
			c = fgetc(file);
		if (c == 'v')
		{
			fscanf(file, "%lf%lf%lf\n", &x, &y, &z);
			vertices[vertex_num].update(vertex_num, x, y, z);
			vertex_num++;
		}
		else if (c == 'f')
		{
			fscanf(file, "%d%d%d\n", &a, &b, &cc);
			a--; b--; cc--;
			faces.push_back(Face(a, b, cc));
			vertices[a].planes.insert(face_num);
			vertices[b].planes.insert(face_num);
			vertices[cc].planes.insert(face_num);
			vertices[a].adjacent_v.insert(b);
			vertices[a].adjacent_v.insert(cc);
			vertices[b].adjacent_v.insert(a);
			vertices[b].adjacent_v.insert(cc);
			vertices[cc].adjacent_v.insert(a);
			vertices[cc].adjacent_v.insert(b);
			face_num++;
		}
		c = 0;
	}

	fclose(file);
	goal_num = face_num * rate;
	start_check = face_num * (0.99 - rate) / 2;
	for (int i = 0;i < vertex_num;i++)
		vertices[i].calQ();
}

void simplify()
{
	heapify();
	clock_t start = clock();
	cout << "start" << endl;
	
	for (int i = 0; ;i += 1)
	{
		if (i % 10000 == 0)
		{
			clock_t end = clock();
			cout << i << " time used: " << (end - start) / CLOCKS_PER_SEC << endl;
		}
		
		if (i > start_check && i % 100 == 0)
		{
			vector<Face>::iterator it;
			set<Face> face;
			for (it = faces.begin();it != faces.end();it++)
			{
				if (it->v[0] != it->v[1] && it->v[1] != it->v[2] && it->v[0] != it->v[2])
					face.insert(*it);
			}
			
			if (face.size() <= goal_num)
				break;
		}
		Edge e = heap.pop();
		Vertex* v1 = vertices + e.v1;
		Vertex* v2 = vertices + e.v2;
		//create new vertex
		Vertex* vnew = vertices + vertex_num;
		vnew->update(vertex_num, e.vbar[0], e.vbar[1], e.vbar[2]);
		vnew->calQ(e.v1, e.v2);
		heap.remove(e);

		//create adjacent vertices set
		set<int> ad_v;
		ad_v.clear();
		set<int>::iterator it;
		for (it = v1->adjacent_v.begin();it != v1->adjacent_v.end();it++)
		{
			if (*it == v2->id) continue;
			heap.remove(*it, v1->id);
			vertices[*it].adjacent_v.erase(v1->id);	
			ad_v.insert(*it);
		}
		for (it = v2->adjacent_v.begin();it != v2->adjacent_v.end();it++)
		{
			if (*it == v1->id) continue;
			heap.remove(*it, v2->id);
			vertices[*it].adjacent_v.erase(v2->id);
			ad_v.insert(*it);
		}
		for (it = ad_v.begin();it != ad_v.end();it++)
		{
			vnew->adjacent_v.insert(*it);
			vertices[*it].adjacent_v.insert(vertex_num);
		}
		v1->valid = false;
		v2->valid = false;
		
		for (it = v1->planes.begin();it != v1->planes.end();it++)
		{
			vnew->planes.insert(*it);
			for (int j = 0;j < 3;j++)
				if (faces[*it].v[j] == v1->id)
					faces[*it].v[j] = vertex_num;
		}
		for (it = v2->planes.begin();it != v2->planes.end();it++)
		{
			vnew->planes.insert(*it);
			for (int j = 0;j < 3;j++)
				if (faces[*it].v[j] == v2->id)
					faces[*it].v[j] = vertex_num;
		}


		//create new edge
		for (it = ad_v.begin();it != ad_v.end();it++)
		{
			Edge e(*it, vertex_num);
			e.update();
			heap.push(e);
		}
		vertex_num++;
	}
}

void heapify()
{
	//put all edges into heap
	for (int i = 0;i < vertex_num;i++)
	{
		Vertex* v = vertices + i;
		set<int>::iterator it;
		for (it = v->adjacent_v.begin();it != v->adjacent_v.end();it++)
		{
			if (i < *it)	break;
			Edge e(*it, i);
			e.update();
			heap.push(e);
		}
	}
}

void output(string filename)
{
	FILE* f;
	/*if (filename != "")
		f = fopen(filename);*/
	f = fopen("D:/models/bunny_0.03.obj", "w");
	
	//output vertices
	int count = 0;
	for (int i = 0;i < vertex_num;i++)
	{
		Vertex* v = vertices + i;
		if (!v->valid) continue;
		count++;
		v->id = count;
		fprintf(f, "v %lf %lf %lf\n", v->coordinate[0], v->coordinate[1], v->coordinate[2]);
	}

	//output faces
	vector<Face>::iterator it;
	set<Face> face;
	for (it = faces.begin();it != faces.end();it++)
	{
		if (it->v[0] != it->v[1] && it->v[1] != it->v[2] && it->v[0] != it->v[2])
			face.insert(*it);
	}
	for (set<Face>::iterator it=face.begin();it!=face.end();it++)
		fprintf(f, "f %d %d %d\n", vertices[it->v[0]].id, vertices[it->v[1]].id, vertices[it->v[2]].id);

}


Matrix::Matrix(double a00, double a01, double a02, double a03, double a10, double a11, double a12, double a13,
	double a20, double a21, double a22, double a23, double a30, double a31, double a32, double a33)
{
	data[0][0] = a00; data[0][1] = a01; data[0][2] = a02; data[0][3] = a03;
	data[1][0] = a10; data[1][1] = a11; data[1][2] = a12; data[1][3] = a13;
	data[2][0] = a20; data[2][1] = a21; data[2][2] = a22; data[2][3] = a23;
	data[3][0] = a30; data[3][1] = a31; data[3][2] = a32; data[3][3] = a33;
}

void Vertex::calQ()
{
	Vertex *v2, *v3;
	//v1 = vertices + v;
	set<int>::iterator it1, it2;
	double a, b, c, d;
	a = b = c = d = 0;
	double a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33;
	a00 = a01 = a02 = a03 = a10 = a11 = a12 = a13 = a20 = a21 = a22 = a23 = a30 = a31 = a32 = a33 = 0;
	//Matrix* Q = new Matrix;
	for (it1 = adjacent_v.begin();it1 != adjacent_v.end();it1++)
	{
		for (it2 = adjacent_v.begin();it2 != adjacent_v.end();it2++)
		{
			if ((*it1) < (*it2) && isEdge(*it1, *it2))
			{
				v2 = vertices + (*it1);
				v3 = vertices + (*it2);
				double x1, x2, y1, y2, z1, z2;
				x1 = coordinate[0] - v2->coordinate[0];
				x2 = coordinate[0] - v3->coordinate[0];
				y1 = coordinate[1] - v2->coordinate[1];
				y2 = coordinate[1] - v3->coordinate[1];
				z1 = coordinate[2] - v2->coordinate[2];
				z2 = coordinate[2] - v3->coordinate[2];
				double A = y1 * z2 - y2 * z1;
				double B = -x1 * z2 + x2 * z1;
				double C = x1 * y2 - x2 * y1;
				double e = 1 / sqrt(A*A + B * B + C * C);
				A *= e;	B *= e; C *= e;
				double D = -(A*coordinate[0] + B * coordinate[1] + C * coordinate[2]);
				a += A;	b += B;	c += C;		d += D;
				a00 += A * A; a01 += A * B; a02 += A * C; a03 += A * D;
				a10 += B * A; a11 += B * B; a12 += B * C; a13 += B * D;
				a20 += C * A; a21 += C * B; a22 += C * C; a23 += C * D;
				a30 += D * A; a31 += D * B; a32 += D * C; a33 += D * D;
			
			}
		}
	}
	Q[0][0] = a00; Q[0][1] = a01; Q[0][2] = a02; Q[0][3] = a03;
	Q[1][0] = a10; Q[1][1] = a11; Q[1][2] = a12; Q[1][3] = a13;
	Q[2][0] = a20; Q[2][1] = a21; Q[2][2] = a22; Q[2][3] = a23;
	Q[3][0] = a30; Q[3][1] = a31; Q[3][2] = a32; Q[3][3] = a33;
}

void Vertex::calQ(int a1, int a2)
{
	Vertex* v1 = vertices + a1;
	Vertex* v2 = vertices + a2;
	for (int i = 0;i < 4;i++)
		for (int j = 0;j < 4;j++)
			Q[i][j] = v1->Q[i][j] + v2->Q[i][j];
}

void Vertex::update(int i, double c1, double c2, double c3)
{
	id = i;
	coordinate[0] = c1;	coordinate[1] = c2;	coordinate[2] = c3;
}

Vertex::Vertex(int i, double c1, double c2, double c3): id(i), valid(true)
{
	coordinate[0] = c1;
	coordinate[1] = c2;
	coordinate[2] = c3;
}

void Edge::update()
{
	Vertex* p1 = vertices + v1;
	Vertex* p2 = vertices + v2;
	Matrix* Q1 = getQ(v1);
	Matrix* Q2 = getQ(v2);
	Matrix* Q = add(Q1, Q2);
	double temp[4];	//save the value of Q[3]	
	for (int i = 0;i < 4;i++)
		temp[i] = Q->data[3][i];

	Q->data[3][0] = Q->data[3][1] = Q->data[3][2] = 0;
	Q->data[3][3] = 1;
	double c = 0;
	Vector_4* v;
	if (abs(det(Q)) > epsilon)		//det(Q) is not zero
	{
		Vector_4 b(0, 0, 0, 1);
		v = solve(Q, &b);
		//restore the value of Q[3]
		for (int i = 0;i < 4;i++)
			Q->data[3][i] = temp[i];
		//calculate cost
		for (int i = 0;i < 4;i++)
		{
			double cc = 0;
			for (int j = 0;j < 4;j++)
				cc += v->data[j] * Q->data[i][j];
			c += cc * v->data[i];
		}
	}
	else  //vbar = (v1+v2)/2
	{
		v = new Vector_4();
		Vector_4 vnew[3];
		for (int i = 0;i < 3;i++)		
		{
			vnew[0].data[i] = (p1->coordinate[i] + p2->coordinate[i]) / 2;		//vnew[0]=(v1+v2)/2
			vnew[1].data[i] = p1->coordinate[i];		//vnew[1]=v1
			vnew[2].data[i] = p2->coordinate[i];		//vnew[2]=v2
			vnew[i].data[3] = 1;
		}
		v->data[3] = 1;
		//restore the value of Q[3]
		for (int i = 0;i < 4;i++)
			Q->data[3][i] = temp[i];
		//calculate cost
		double cost1, cost2, cost3;
		cost1 = cost2 = cost3 = 0;
		for (int i = 0;i < 4;i++)
		{
			double cc = 0;
			for (int j = 0;j < 4;j++)
				cc += vnew[0].data[j] * Q->data[i][j];
			cost1 += cc * vnew[0].data[i];
		}
		for (int i = 0;i < 4;i++)
		{
			double cc = 0;
			for (int j = 0;j < 4;j++)
				cc += vnew[1].data[j] * Q->data[i][j];
			cost2 += cc * vnew[1].data[i];
		}
		for (int i = 0;i < 4;i++)
		{
			double cc = 0;
			for (int j = 0;j < 4;j++)
				cc += vnew[2].data[j] * Q->data[i][j];
			cost3 += cc * vnew[2].data[i];
		}
		double temp = cost1 < cost2 ? cost1 : cost2;
		c = temp < cost3 ? temp : cost3;
		if (c == cost1)
		{
			for (int i = 0;i < 3;i++)
				v->data[i] = vnew[0].data[i];
		}
		else if (c == cost2)
		{
			for (int i = 0;i < 3;i++)
				v->data[i] = vnew[1].data[i];
		}
		else
		{
			for (int i = 0;i < 3;i++)
				v->data[i] = vnew[2].data[i];
		}
	}	
	cost = c;

	for (int i = 0;i < 3;i++)
		vbar[i] = v->data[i];
	delete Q;
	delete Q1;
	delete Q2;
	delete v;
}

void Heap::push(Edge & e)
{
	e.id = size;
	int a1 = e.v1 > e.v2 ? e.v1 : e.v2;
	int a2 = e.v1 + e.v2 - a1;
	isDeleted[size] = false;
	ver2edge[make_pair(a1, a2)] = size++;
	queue.push(e);
}

void Heap::remove(Edge e)
{
	int a1 = e.v1 > e.v2 ? e.v1 : e.v2;
	int a2 = e.v1 + e.v2 - a1;
	int id = ver2edge[make_pair(a1, a2)];
	isDeleted[id] = true;
}

void Heap::remove(int a, int b)
{
	int v1 = a > b ? a : b;
	int v2 = a + b - v1;
	if (ver2edge.count(make_pair(v1, v2)))
		isDeleted[ver2edge[make_pair(v1, v2)]] = true;
}

Edge Heap::pop()
{
	while (isDeleted[queue.top().id] || !vertices[queue.top().v1].valid || !vertices[queue.top().v2].valid)
		queue.pop();
	Edge e = queue.top();
	queue.pop();
	return e;
}

bool Face::operator<(const Face & f) const
{
	if (v[0] < f.v[0])
		return true;
	else if (v[0] > f.v[0])
		return false;
	if (v[1] < f.v[1])
		return true;
	else if (v[1] > f.v[1])
		return false;
	return v[2] < f.v[2];
}
