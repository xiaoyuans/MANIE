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

/*
 * �ó�����01�������ѹ����֪ʹ��
 * ����ʾ��ͼ�ο� ���ס�Reconstructing Propagation Networks with Natural
 * Diversity and Identifying Hidden Source��
 *
 * �ó���ɽ���4��������
 * '-s' �������ʼ�ڵ㣬 Ĭ����0  -s 10
 * '-t' ����ڵ㷶Χ  Ĭ����ȫ�� -s 10 -t 10 ָ����10~19�Žڵ�
 * '-s1' ѡbasis��hamming distance ���ƶȣ�ѡ��ͬ�����ƶ�
 * '-s2' ��ÿ��base ����ͬ�����ƶ�
*/

#ifdef _WIN32
//#define strcasecmp _stricmp
//#define strncasecmp _strnicmp
#endif
#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#endif

using namespace std;
typedef vector<vector<int> > Mat;
Mat input();//���ڶ�ȡ�����ļ����ļ���ʽΪ01����
int main(int argc, char** argv)
{
	//���������������������
	int numt, flag, inumt, ham, jnumt, Aibs, asize, yy, lidof0 = 0, cbn = 10000, nstart = 0;
	char ama[60] = "Amatrix.txt", yma[60] = "Ymatrix.txt";
	char* rx;
	Mat aa = input();
	numt = aa.size();
	const int totalN = aa[0].size();
	int nend = totalN;
	float sim1 = round(0.25 * totalN);
	float sim2 = round(0.45 * totalN);

	//���������ָ������������ĸ��ڵ������
	for (int i = 1; i < argc; i++)
	{
		if (!strcasecmp(argv[i], "-s"))
		{
			nstart = atoi(argv[i + 1]);
			rx = argv[i + 1];
			strcat(ama, rx);
			strcat(yma, rx);
		}
		else if (!strcasecmp(argv[i], "-t"))
		{
			nend = atoi(argv[i + 1]);
			nend += nstart;
		}
		else if (!strcasecmp(argv[i], "-s1"))
		{
			sim1 = atof(argv[i + 1]);
			sim1 = sim1 * totalN;
		}
		else if (!strcasecmp(argv[i], "-s2"))
		{
			sim2 = atof(argv[i + 1]);
			sim2 = sim2 * totalN;
		}
	}
	if (nstart >= totalN)
	{
		return 0;
	}
	else
	{
		int* Abasis;//�洢basis
		int* Aidof0;//�洢0��λ��
		int* XY;
		int* Aidofbs;
		int* Aidsbs;
		// from (1,-1) to (1,0)
		bool** a;//�洢01����
		a = (bool**)malloc(numt * sizeof(bool*));
		Aidof0 = (int*)malloc(numt * sizeof(int));
		Aidofbs = (int*)malloc(numt * sizeof(int));
		Aidsbs = (int*)malloc(numt * sizeof(int));
		Abasis = (int*)malloc(totalN * sizeof(int));
		XY = (int*)malloc(totalN * sizeof(int));
		for (int i = 0; i < numt; i++)
		{
			a[i] = (bool*)malloc(totalN * sizeof(bool));
			for (int j = 0; j < totalN; j++)
			{
				//~ a[i][j] = (abs(aa[i][j])+aa[i][j])/2;
				a[i][j] = aa[i][j];
			}
		}
		aa.clear();
		time_t time0, time1;
		time0 = time(&time0);

		cout << sim1 << ' ' << sim2 << endl;
		if (nend > totalN) nend = totalN;
		ofstream onf1(ama);//����ļ� \Phi
		ofstream onf2(yma);//����ļ� Y

		//��������
		for (int node = nstart; node < nend; node++)
		{
			cout << "node: " << node << flush;
			// find index of si = 0
			// Ѱ��01,00���洢0��λ�ã�����Aidof0
			lidof0 = 0;
			for (int i = 0; i < (numt - 1); i++)
			{
				if (a[i][node] == 0)
				{
					Aidof0[lidof0] = i;
					lidof0++;
				}
			}
			/*           find basis             */
			//��һ�����ҵ����ʵ�basis
			cbn = lidof0;
			for (int i = 0; i < cbn; i++)
			{
				Aidofbs[i] = Aidof0[i];
			}
			//~ float xx1 = 0.35;
			flag = totalN;
			jnumt = cbn;
			while (flag) {
				inumt = jnumt;

				if (inumt == 0)
				{
					cout << "flag: " << flag << endl << flush;
					int iflag = totalN - flag;
					while (flag)
					{
						Abasis[totalN - flag] = Abasis[flag % iflag];
						--flag;
					}
					break;
				}

				// calculate the candidate base,leave the useful base
				// ����ɹ�ѡ���basis
				jnumt = 0;
				for (int i = 1; i < inumt; i++)
				{
					ham = 0;
					Aibs = Aidofbs[0];
					for (int j = 0; j < totalN; j++)
					{
						ham += a[Aidofbs[i]][j] ^ a[Aibs][j];
					}
					if (ham > sim1)//leave the useful base�����º���������basis
					{
						Aidsbs[jnumt] = Aidofbs[i];
						jnumt++;
					}
				}
				asize = jnumt;
				for (int i = 0; i < asize; i++)
				{
					Aidofbs[i] = Aidsbs[i];
				}
				Abasis[totalN - flag] = Aidofbs[0];
				--flag;
			}
			cout << "  flag: " << flag << endl;
			//~ time1 = time(&time1);
			//~ cout<<difftime(time0,time1)<<endl;


			//�ڶ�����������ÿ��base���Ƶ����У�ͳ��ƽ��ֵ�����������
			inumt = lidof0;
			for (int ibs = 0; ibs < totalN; ibs++) {
				//~ cout<<ibs<<" "<<flush;
				//~ xx1 = 0.45;
				jnumt = 0;
				yy = 0;
				for (int i = 0; i < totalN; i++) {
					XY[i] = 0;
				}
				for (int i = 0; i < inumt; i++) {
					ham = 0;
					for (int j = 0; j < totalN; ++j) {
						ham += a[Aidof0[i]][j] ^ a[Abasis[ibs]][j];

					}
					if (ham < sim2)//ͳ�ƺ��������С�����е�ƽ��ֵ
					{
						//Aidsbs[jnumt]=Aidof0[i];
						yy += a[Aidof0[i] + 1][node];
						for (int j = 0; j < totalN; ++j) {
							XY[j] += a[Aidof0[i]][j];
						}
						jnumt++;
					}
				}
				//~ number = jnumt;
				for (int j = 0; j < totalN; ++j)
				{
					onf1 << (1.0 * XY[j] / jnumt) << " ";//���\Phi
				}
				//~ onf2<<log(1.0-1.0*yy/jnumt)/log(1-0.2)<<" ";//���Y
				onf2 << -1.0 * log(1.0 - 1.0 * yy / jnumt) << " ";
				onf1 << endl;
			}
			onf2 << endl;
			//~ time1 = time(&time1);
			//~ cout<<difftime(time0,time1)<<endl;
		}
		onf1.close();
		onf2.close();
		//�ͷ��ڴ�
		for (int t = 0; t < numt; t++) {
			free(a[t]);
		}
		free(a);
		free(Abasis);
		free(XY);
		free(Aidof0);
		free(Aidofbs);
		free(Aidsbs);


		time1 = time(&time1);
		cout << difftime(time0, time1);
		return 0;
	}
}



Mat input()
{
	/*
	 *�ú������ڶ���01���ݣ���άvector���ݣ�T��ʱ�䣩��N���ڵ㣩��
	*/
	ifstream in("s_time_state.txt");
	Mat a;
	istringstream istr;
	string str;
	vector<int> tmpvec;
	while (getline(in, str))
	{
		istr.str(str);
		int tmp;
		while (istr >> tmp)
		{
			tmpvec.push_back(tmp);
		}
		a.push_back(tmpvec);
		tmpvec.clear();
		istr.clear();
	}
	cout << a.size() << endl;
	return a;
}
