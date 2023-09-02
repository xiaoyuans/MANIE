//�ó�������ģ���������
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

int N = 115; // number of nodes
int T = 115 * 100; // length of timeseries��run x round add the initial round T = x+1
#define noise  0 // �����ʱ�����еĽڵ������ fraction of unknown nodes
#define trans  0.1 // ״ֵ̬��ת�ĸ���
#define stay 0 //״ֵ̬һֱû�仯�ĸ���//�������֮��Ymatrix�������
int** state; // �ڵ�״̬
int** neighbor;
FILE* fp4, * fp5, * fp6, * fp7; // ����ļ�

void ini_pointer() // �����ڴ�
{
	state = (int**)malloc(T * sizeof(int));
	for (int i = 0; i < T; i++)
	{
		state[i] = (int*)malloc(N * sizeof(int));
	}
	neighbor = (int**)malloc(N * sizeof(int));
	for (int i = 0; i < N; i++)
	{
		neighbor[i] = (int*)malloc(N * sizeof(int));
	}
}

void free_pointer() // �ͷ��ڴ�
{
	for (int i = 0; i < T; i++)
	{
		free(state[i]);
	}
	free(state);
	for (int i = 0; i < N; i++)
	{
		free(neighbor[i]);
	}
	free(neighbor);
}


int main(int argc, char** argv)
{
	ini_pointer();
	fp4 = fopen("s_time_state.txt", "r");
	if (fp4 == NULL)
	{
		printf("�ļ���ȡ����...");
		return -1;
	}
	for (int i = 0; i < T; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fscanf(fp4, "%d", &state[i][j]);/*ÿ�ζ�ȡһ������fscanf���������ո���߻��н���*/
		}
	}
	fclose(fp4);

	int min = 0;
	int max = N - 1;
	int num = (N * noise);
	int num_1 = (N * stay);
	int rnd;
	int flag;
	double temp;
	vector<int> diff;
	vector<int> diff_1;
	vector<int> tmp;//�洢ʣ�����
	//��ʼ��
	for (int i = min; i < max + 1; i++)
	{
		tmp.push_back(i);
	}
	srand((unsigned)time(0)); //��ʼ�����������
	for (int i = 0; i < num; i++)
	{
		do {
			rnd = min + rand() % (max - min + 1);

		} while (tmp.at(rnd - min) == -1);
		diff.push_back(rnd);
		tmp.at(rnd - min) = -1;
	}

	for (int i = min; i < max + 1; i++)
	{
		tmp.push_back(i);
	}
	srand((unsigned)time(0)); //��ʼ�����������
	for (int i = 0; i < num_1; i++)
	{
		do {
			rnd = min + rand() % (max - min + 1);

		} while (tmp.at(rnd - min) == -1);
		diff_1.push_back(rnd);
		tmp.at(rnd - min) = -1;
	}

	// cout << diff[0] << endl;
	//cout << state[0][1] << endl;

	for (int i = 0; i < T; i++)
	{
		for (int j = 0; j < N; j++)
		{
			temp = (double)rand() / RAND_MAX;
			if (temp < trans)
			{
				state[i][j] = 1 - state[i][j];
			}
		}
		for (int k = 0; k < num; k++)
		{
			state[i][diff[k]] = 0;
		}
		for (int k = 0; k < num_1; k++)
		{
			state[i][diff_1[k]] = state[0][diff_1[k]];
		}
	}

	fp5 = fopen("s_time_state_noise.txt", "w");
	for (int i = 0; i < T; i++)
	{
		for (int j = 0; j < N; j++)
		{
			flag = 1;
			for (int k = 0; k < num; k++)
			{
				if (j == diff[k])
				{
					flag = 0;
					break;
				}

			}
			if (flag)
			{
				fprintf(fp5, "%d ", state[i][j]);
			}		
		}
		fprintf(fp5, "\n");
	}
	fclose(fp5);

	fp6 = fopen("Xmatrix.txt", "r");
	if (fp6 == NULL)
	{
		printf("�ļ���ȡ����...");
		return -1;
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fscanf(fp6, "%d", &neighbor[i][j]);/*ÿ�ζ�ȡһ������fscanf���������ո���߻��н���*/
		}
	}
	fclose(fp6);

	fp7 = fopen("Xmatrix_noise.txt", "w");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			flag = 1;
			for (int k = 0; k < num; k++)
			{
				if ((i == diff[k]) ||  (j == diff[k]))
				{
					flag = 0;
					break;
				}

			}
			if (flag)
			{
				fprintf(fp7, "%d ", neighbor[i][j]);
			}
		}
		fprintf(fp7, "\n");
	}
	fclose(fp7);

	free_pointer();
	cout << "Done!";
	return 0;
}

/*
int main()
{
	vector<string> vv;
	ifstream fin("c_time_state.txt");
	if (!fin.is_open())
	{
		cout << "open error!" << endl;
	}
	//�����ݴ���vv����(���ַ�����ʽ)
	string temp;
	while (getline(fin, temp))
	{
		vv.push_back(temp);
	}
	return 0;
	*/
 /*
vector<int> GenerateDiffNumber(int min, int max, int num)
{
	int rnd;
	vector<int> diff;
	vector<int> tmp;//�洢ʣ�����
	//��ʼ��
	for (int i = min; i < max + 1; i++)
	{
		tmp.push_back(i);
	}
	srand((unsigned)time(0)); //��ʼ�����������
	for (int i = 0; i < num; i++)
	{
		do {
			rnd = min + rand() % (max - min + 1);

		} while (tmp.at(rnd - min) == -1);
		diff.push_back(rnd);
		tmp.at(rnd - min) = -1;
	}
	return diff;
}
*/