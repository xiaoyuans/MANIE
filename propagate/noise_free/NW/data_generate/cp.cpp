/*
 * 该程序用于输出cp过程时间序列
 * 可接收的参数
 * -N 节点数量
 * -T 时间序列长度
 * -D 节点是否异质
*/

#include <iostream>
#include <vector>
#include <deque>
#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <sstream>
using namespace std;

#define iniP 0.5  //initial partition of I individuals
#define infP 0.7    // probability of getting virus
#define recP 0.2  // probability of recover from I

#ifdef _WIN32
//#define strcasecmp _stricmp
//#define strncasecmp _strnicmp
#endif
#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#endif

int N = 40; // number of nodes
int T = 50000; // length of timeseries，run x round add the initial round T = x+1
int diver = 1; // 1 individual has different infect(recover) rate, 0 homogeneous
float noise = 0.0; // 不输出时间序列的节点比例， fraction of unknown nodes
int** nbhmatrix; //存储邻居关系
int** adjmatrix; //存储邻接矩阵
int* known; // 存储是否输入该节点状态的列表，1 表示输出， 0 不输出
int* degree;// 存储每个节点的度， degree of each node
int* nowstate; // 节点当前状态
int* tempstate; // 节点状态临时存储
float* lambda; // 节点染病概率  infect rate of each node
float* delta; // 节点恢复概率  recover rate of each node
FILE* fp, * fp1, * fp2, * fp3; // 输出文件


void ini_pointer() // 分配内存
{
	adjmatrix = (int**)malloc(N * sizeof(int*));
	nbhmatrix = (int**)malloc(N * sizeof(int*));
	for (int i = 0; i < N; i++)
	{
		adjmatrix[i] = (int*)malloc(N * sizeof(int));
		nbhmatrix[i] = (int*)malloc(N * sizeof(int));
	}
	nowstate = (int*)malloc(N * sizeof(int));
	tempstate = (int*)malloc(N * sizeof(int));
	degree = (int*)malloc(N * sizeof(int));
	lambda = (float*)malloc(N * sizeof(float));
	delta = (float*)malloc(N * sizeof(float));
}
void free_pointer() // 释放内存
{
	free(nowstate);
	free(tempstate);
	free(degree);
	free(lambda);
	for (int i = 0; i < N; i++)
	{
		free(adjmatrix[i]);
		free(nbhmatrix[i]);
	}
	free(adjmatrix);
	free(nbhmatrix);
	free(delta);
}
void initialize() // 初始化各指针变量
{
	int i;
	double temp;
	FILE* bp;
	FILE* lp;
	bp = fopen("lambda.txt", "w");
	lp = fopen("delta.txt", "w");
	for (i = 0; i < N; i++)
	{
		tempstate[i] = 0.0;

		if (diver == 1)
		{
			temp = (double)rand() / RAND_MAX;
			lambda[i] = infP + temp * 0.2;
			temp = (double)rand() / RAND_MAX;
			delta[i] = recP + temp * 0.2;
		}
		else
		{
			lambda[i] = infP + 0.1;
			delta[i] = recP + 0.1;
		}
		fprintf(bp, "%f ", lambda[i]);
		fprintf(lp, "%f ", delta[i]);
		temp = (double)rand() / RAND_MAX;
		if (temp < iniP)
		{
			nowstate[i] = 1;
		}
		else
		{
			nowstate[i] = 0;
		}
	}
	fclose(bp);
	fclose(lp);

}

void input_adj() // 读入节点邻接关系
{
	int m, n;
	int i, j;
	for (i = 0; i < N; i++)
	{
		degree[i] = 0;
	}
	fp = fopen("adj_nw.txt", "r");
	while (fscanf(fp, "%d %d", &m, &n) == 2)
	{
		printf("%d %d\n", m, n);
		nbhmatrix[m][degree[m]] = n;
		degree[m] ++;
	}

	fclose(fp);
	fp3 = fopen("Xmatrix.txt", "w");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			adjmatrix[i][j] = 0;
		}
	}
	int neighbor;
	for (i = 0; i < N; i++) {
		for (j = 0; j < degree[i]; j++) {
			neighbor = nbhmatrix[i][j];
			adjmatrix[i][neighbor] = 1;
		}
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			fprintf(fp3, "%d ", adjmatrix[i][j]);   // N by N network nbhmatrix matrix
		}
		fprintf(fp3, "\n");
	}
	fclose(fp3);
}


int communication(int timeT)
{
	int i, j;
	double temp, tload = 0.0; // tload is the temp receiving load from a neighbor
	int neighbor;
	for (i = 0; i < N; i++)
	{ // node loop  should add neighbors loop
		// just choose one neighbor
		tload = 0.0;
		if (nowstate[i] < 1.0)
		{//for S node may get virus

			for (j = 0; j < degree[i]; j++)
			{
				neighbor = nbhmatrix[i][j];
				tload += nowstate[neighbor];
			}
			/*
			j=rand()%degree[i];
			neighbor = nbhmatrix[i][j];
			tload = nowstate[neighbor];
			*/
			temp = (float)rand() / RAND_MAX;
			if (temp < lambda[i] * tload / (float)degree[i])
				//if(temp < P*tload)
			{ // get virus from neighbors
				tempstate[i] = 1.0;
			}
			else
			{
				tempstate[i] = 0.0;
			}
		}
		else
		{//for I node won't get virus but recover
			temp = (float)rand() / RAND_MAX;
			if (temp < delta[i])
			{// I individuals have a chance to recovery
				tempstate[i] = 0.0;
			}
			else
			{
				tempstate[i] = 1.0;
			}
		}
	} // node loop over
	float sum = 0.0;
	for (i = 0; i < N; i++)
	{
		nowstate[i] = tempstate[i];
		sum += nowstate[i];
	}
	int Node = N;
	sum = sum / (float)Node;
	//tsum=(tsum*(float)(timeT-1)+sum)/(float)timeT;
	//printf("%f,%f\n",tsum,sum);
	/*
	FILE *f;
	f=fopen("tt.txt","a");
	fprintf(f,"%f\n",sum);
	fclose(f);
	*/
	if (sum > 0.00001)
	{
		return 1;
	}
	else
	{
		return 0;// 表示当前轮所有的节点状态都为0
	}
}
void output()
{//output the timeseries state of each node
	int i;
	for (i = 0; i < N; i++)
	{// node loop
		fprintf(fp, "%d ", (int)nowstate[i]);
	}
	fprintf(fp, "\n");
}


int main(int argc, char** argv)
{
	time_t   t;
	srand((unsigned)time(&t));
	int timeT = 0, flag = 1;
	//argc 表示传入main函数的参数个数
	//argv 表示传入main函数的参数序列或指针，并且第一个参数argv[0]一定是程序的名称，并且包含了程序所在的完整路径
	// 所以确切的说需要我们输入的main函数的参数个数应该是argc-1个；
	for (int i = 1; i < argc; i++)
	{
		if (!strcasecmp(argv[i], "-T"))
		{
			T = atoi(argv[i + 1]);//atoi()函数将数字格式的字符串转换为整数类型。
		}
		else if (!strcasecmp(argv[i], "-N"))
		{
			N = atoi(argv[i + 1]);
		}
		else if (!strcasecmp(argv[i], "-D"))
		{
			diver = atoi(argv[i + 1]);
		}
	}

	ini_pointer();
	input_adj();
	initialize();
	fp = fopen("c_time_state.txt", "w"); // write A matrices
	while (flag && (timeT < T))
	{
		output();
		flag = communication(timeT);
		timeT++;
	}
	fclose(fp);
	free_pointer();
	cout << timeT << endl;
	cout << flag << "hello world!";
	return 0;
}
