/*
 * �ó����������sis����ʱ������
 * �ɽ��յĲ���
 * -N �ڵ�����
 * -T ʱ�����г���
 * -D �ڵ��Ƿ�����
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

#ifdef _WIN32
//#define strcasecmp _stricmp
//#define strncasecmp _strnicmp
#endif
#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#endif

#define iniP 0.5  //initial partition of I individuals
#define infP 0.2    // probability of getting virus
#define recP 0.4  // probability of recover from I

int N = 34; // number of nodes
int T = 50000; // length of timeseries��run x round add the initial round T = x+1
int diver = 1; // 1 individual has different infect(recover) rate, 0 homogeneous
float noise = 0.0; // �����ʱ�����еĽڵ������ fraction of unknown nodes
int** nbhmatrix; //�洢�ھӹ�ϵ
int** adjmatrix; //�洢�ڽӾ���
int* known; // �洢�Ƿ�����ýڵ�״̬���б���1 ��ʾ����� 0 �����
int* degree;// �洢ÿ���ڵ�Ķȣ� degree of each node
int* nowstate; // �ڵ㵱ǰ״̬
int* tempstate; // �ڵ�״̬��ʱ�洢
float* lambda; // �ڵ�Ⱦ������  infect rate of each node
float* delta; // �ڵ�ָ�����  recover rate of each node
FILE* fp, * fp1, * fp2, * fp3; // ����ļ�


void ini_pointer() // �����ڴ�
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
void free_pointer() // �ͷ��ڴ�
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
void initialize() // ��ʼ����ָ�����
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

void input_adj() // ����ڵ��ڽӹ�ϵ
{
	int m, n;
	int i, j;
	for (i = 0; i < N; i++)
	{
		degree[i] = 0;
	}
	fp = fopen("adj.txt", "r");
	while (fscanf(fp, "%d %d", &m, &n) == 2)
	{
		//printf("%d %d\n",m,n);
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


int communication(int timeT)// �ڵ�״̬��Ⱦ���̣� ���״̬��ȫΪ0�򷵻�1
{
	int i, j, tload, sum;
	double temp; // tload is the temp receiving load from a neighbor
	int neighbor;

	for (i = 0; i < N; i++)
	{ // for each node tranverse its neighbors

		tload = 0;
		if (nowstate[i] < 1)
		{//for S node may get virus
			for (j = 0; j < degree[i]; j++)
			{
				temp = (double)rand() / RAND_MAX;
				if (temp < lambda[i])
				{ // get virus from neighbors
					neighbor = nbhmatrix[i][j];
					tload += nowstate[neighbor];//tload should be state of S or I  marked as 0 and 1,then Range shows its own state now,thus Range = nowstate[i][t-1]
				}
			}
			//nowstate[i][t]=min(1,tempstate[i][t];
			if (tload > 0)
			{
				tempstate[i] = 1;
			}
			else
			{
				tempstate[i] = 0;
			}
		}
		else
		{//for I node won't get virus but recover
			temp = (double)rand() / RAND_MAX;
			if (temp < delta[i])
			{// I individuals have a chance to recovery
				tempstate[i] = 0;
			}
			else
			{
				tempstate[i] = 1;
			}
		}
	} // node loop over
	sum = 0;
	for (i = 0; i < N; i++)
	{
		nowstate[i] = tempstate[i];
		sum += nowstate[i];
	}
	int Node = N;
	temp = (float)sum / (float)Node;
	//printf("%f\n",temp);
	if (temp < 0.0001)
	{
		return 0;
	}
	else
	{
		return 1;
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

	for (int i = 1; i < argc; i++)
	{
		if (!strcasecmp(argv[i], "-T"))
		{
			T = atoi(argv[i + 1]);
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
	fp = fopen("s_time_state.txt", "w"); // write A matrices
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