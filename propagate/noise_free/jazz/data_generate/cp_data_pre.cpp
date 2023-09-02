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
 * 该程序处理01串供耦合压缩感知使用
 * 处理示意图参考 文献‘Reconstructing Propagation Networks with Natural
 * Diversity and Identifying Hidden Source’
 *
 * 该程序可接受2个外界参数
 * '-s' 计算的起始节点， 默认是0  -s 10
 * '-t' 计算节点范围  默认是全部 -s 10 -t 10 指计算10~19号节点
 * '-s1' 选basis的hamming distance 相似度，选不同的相似度
 * '-s2' 给每个base 找相同的相似度
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
Mat input();//用于读取数据文件，文件格式为01数组
int main(int argc, char** argv)
{
	//定义各变量，并读入数据
	int numt, flag, inumt, ham, jnumt, Aibs, asize, yy, lidof0 = 0, cbn = 10000, nstart = 0;
	char ama[60] = "Amatrix.txt", yma[60] = "Ymatrix.txt";
	char* rx;
	Mat aa = input();
	numt = aa.size();
	const int totalN = aa[0].size();
	int nend = totalN;
	float sim1 = round(0.35 * totalN);
	float sim2 = round(0.45 * totalN);

	cout << totalN << endl;
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
		int* Abasis;//存储basis
		int* Aidof0;//存储0的位置
		int* XY;
		int* Aidofbs;
		int* Aidsbs;
		// from (1,-1) to (1,0)
		bool** a;//存储01数组
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
		//读入参数，指定本程序计算哪个节点的序列

		cout << sim1 << ' ' << sim2 << endl;
		if (nend > totalN) nend = totalN;
		ofstream onf1(ama);//输出文件 \Phi
		ofstream onf2(yma);//输出文件 Y

		//处理数据
		for (int node = nstart; node < nend; node++)
		{
			cout << "node: " << node << flush;
			// find index of si = 0
			// 寻找01,00，存储0的位置，存入Aidof0
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
			//第一步，找到合适的basis
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
				// 计算可供选择的basis
				jnumt = 0;
				for (int i = 1; i < inumt; i++)
				{
					ham = 0;
					Aibs = Aidofbs[0];
					for (int j = 0; j < totalN; j++)
					{
						ham += a[Aidofbs[i]][j] ^ a[Aibs][j];
					}
					if (ham > sim1)//leave the useful base，留下海明距离大的basis
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


			//第二步，计算与每个base相似的序列，统计平均值，并输出数据
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
					if (ham < sim2)//统计海明距离较小的序列的平均值
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
					onf1 << (1.0 * XY[j] / jnumt) << " ";//输出\Phi
				}
				onf2 << (1.0 * yy / jnumt) << " ";//输出Y
				onf1 << endl;
			}
			onf2 << endl;
			//~ time1 = time(&time1);
			//~ cout<<difftime(time0,time1)<<endl;
		}
		onf1.close();
		onf2.close();
		//释放内存
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
	 *该函数用于读入01数据，二维vector数据，T（时间）行N（节点）列
	*/
	ifstream in("c_time_state.txt");
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
