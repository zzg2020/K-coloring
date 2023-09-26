

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include<math.h>
#include<conio.h>
#include<ctype.h>
#define VN 2000
int t_limit;
int seed;
//#define SHOWDEBUG

using namespace std;

string File_Name;
FILE* fp;
char outfilename[30];
int Seed;
int Total_Iterations;
int Max_Non_Improv_One = 2000000;
double starting_time, finishing_time, once_total_time, total_time;

typedef struct Adjacent_Matrix {
    int neighbor;
    struct Adjacent_Matrix* next;
}Adjacent;
Adjacent** A_Matrix;


typedef struct NeighborStep {
    int v; // vertex number
    int sc; // src color
    int dc; // dst color

    NeighborStep(int nv, int srccolor, int dstcolor)
    {
        v = nv;
        sc = srccolor;
        dc = dstcolor;
    }
    NeighborStep()
    {
        v = -1;
        sc = -1;
        dc = -1;
    }
}Nstep;

int N, K, G_K;  // node number and color number
int f, f_best, global_f_best;
int* Color; // color array for each vertex
int* Best_Color, * Global_Best_Color; // color array for each vertex
int** Adjcaent_Color_Table;  // incremental matrix for each vertex being colored each color
int** Edge;   // adjacent matrix
int** TabuTenure;  // tabu tenure



/*****************************************************************************/
/*****************          1. Initializing           ************************/
/*****************************************************************************/
//1.1 Intializations
void Initializing()
{
    int i, j, x, y, x1, x2;
    //int nb_vtx=0 ;
    int nb_edg = -1, max_edg = 0;


    cin >> N >> nb_edg >> K;
    Color = new int[N];//存放每个顶点所在的颜色集合
    Best_Color = new int[N];
    Global_Best_Color = new int[N];
    //   Move_Freq = new int[N];

    Edge = new int* [N];
    for (x = 0; x < N; x++)
        Edge[x] = new int[N];

    A_Matrix = new Adjacent * [N];//A_Matrix[i] 指向第i个Vertex对应的所有相邻节点
    for (i = 0; i < N; i++)//矩阵初始化
    {
        A_Matrix[i] = new Adjacent;
        A_Matrix[i]->neighbor = 0;
        A_Matrix[i]->next = NULL;
    }

    Adjcaent_Color_Table = new int* [N];//记录每个顶点相邻的点中每个颜色分别有多少个
    for (x = 0; x < N; x++) Adjcaent_Color_Table[x] = new int[K];//每个顶点都维护着k个数据
    TabuTenure = new int* [N];//禁忌表

    for (x = 0; x < N; x++) TabuTenure[x] = new int[K];

    for (x = 0; x < N; x++)
        for (y = 0; y < N; y++)
        {
            Edge[x][y] = 0;
        }

    for(i=0;i<nb_edg;i++)
    {
        cin >> x1 >> x2;
        // cout << x1 << x2 << endl;
        if (x1 < 0 || x2 < 0 || x1 >= N || x2 >= N)
        {
            exit(0);
        }
        Edge[x1][x2] = Edge[x2][x1] = 1;
        max_edg++;
        // add x2 to x1's neighbor list
        Adjacent* p1, * q1;
        p1 = A_Matrix[x1];
        A_Matrix[x1]->neighbor++;//更新x1的邻居信息
        while (p1->next != NULL)
            p1 = p1->next;
        q1 = new Adjacent;
        q1->neighbor = x2;
        q1->next = NULL;
        p1->next = q1;

        // add x1 to x2's neighbor list
        p1 = A_Matrix[x2];
        A_Matrix[x2]->neighbor++;//更新x2的邻居信息
        while (p1->next != NULL)
            p1 = p1->next;
        q1 = new Adjacent;
        q1->neighbor = x1;
        q1->next = NULL;
        p1->next = q1;
    }
    if (0 && max_edg != nb_edg)
    {
        exit(0);
    }
}



void InitTenure()
{
    for (int x = 0; x < N; x++)
    {
        for (int y = 0; y < K; y++)
        {
            TabuTenure[x][y] = 0;
        }
    }
}

void InitColor()
{
    for (int i = 0; i < N; ++i)
    {
        srand(time(NULL));
        Color[i] = rand() % K;
    }
}

void InitAdjcaent_Color_Table()
{

    Adjacent* p = NULL;

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < K; ++j)
            Adjcaent_Color_Table[i][j] = 0;

    for (int i = 0; i < N; ++i)
    {
        p = A_Matrix[i];
        while (p->next)
        {
            p = p->next;
            Adjcaent_Color_Table[i][Color[p->neighbor]]++;
        }
    }

#ifdef SHOWDEBUG
    cout << "Adjcaent_Color_Table" << endl;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < K; ++j)
            cout << Adjcaent_Color_Table[i][j] << " ";
        cout << endl;
    }
#endif
}

void InitConflictNumber()
{
    f = 0;
    for (int i = 0; i < N; ++i)
    {
        f = f + Adjcaent_Color_Table[i][Color[i]];
    }
    f = f / 2;
    f_best = f;
    global_f_best = f;

#ifdef SHOWDEBUG
    cout << "f: " << f << endl;
#endif
}


int NeighborMotion(Nstep& step, int currentcount)
{
    int Delta = 0, Delta_min = 65535, Delta_T = 0, Delta_T_min = 65535;
    vector<Nstep> steps, steps_T;
    for (int i = 0; i < N; ++i)
    {
        if (Adjcaent_Color_Table[i][Color[i]] != 0)
        {
            int m_color = 0;
            while (m_color < K)
            {
                if (m_color == Color[i])
                {
                    ++m_color;
                    continue;
                }
                if (TabuTenure[i][m_color] > currentcount)  //for tabu NeighborMotion
                {
                    Delta_T = Adjcaent_Color_Table[i][m_color] - Adjcaent_Color_Table[i][Color[i]];
                    if (Delta_T < Delta_T_min)
                    {
                        Delta_T_min = Delta_T;
                        steps_T.clear();
                        steps_T.push_back(Nstep(i, Color[i], m_color));//remember moves(u,i,j), u-th vertex change color-i to color-j;
                    }
                    else if (Delta_T == Delta_T_min)
                    {
                        steps_T.push_back(Nstep(i, Color[i], m_color));//add to result
                    }
                }
                else    //for normal NeighborMotion
                {
                    Delta = Adjcaent_Color_Table[i][m_color] - Adjcaent_Color_Table[i][Color[i]];
                    if (Delta < Delta_min)
                    {
                        Delta_min = Delta;
                        steps.clear();
                        steps.push_back(Nstep(i, Color[i], m_color));//remember moves(u,i,j), u-th vertex change color-i to color-j;
                    }
                    else if (Delta == Delta_min)
                    {
                        steps.push_back(Nstep(i, Color[i], m_color));//add to result
                    }
                }
                ++m_color;
            }
        }
    }

    if (f + Delta_T_min < f_best && Delta_T_min < Delta_min)
    {
        int k = rand();
        k = k % steps_T.size();
        step = steps_T[k];//choose to random one from tabu move
        f = f + Delta_T_min;
        return Delta_T_min;
    }
    else
    {
        if (steps.size() == 0)
        {
            step = Nstep(-1, -1, -1);
        }
        else
        {
            int k = rand();
            k = k % steps.size();
            step = steps[k];//choose to random one from normal move
            f = f + Delta_min;
        }
        return Delta_min;
    }
}

void Update(Nstep& step)
{
    //Update f_best
    if (f < f_best)
    {
        f_best = f;
        for (int i = 0; i < N; ++i)
        {
            Best_Color[i] = Color[i];
        }
    }
    //Update Adjcaent_Color_Table
    Adjacent* p = A_Matrix[step.v];
    while (p->next)
    {
        p = p->next;
        Adjcaent_Color_Table[p->neighbor][step.sc] -= 1;
        Adjcaent_Color_Table[p->neighbor][step.dc] += 1;
    }
}

bool GetParam(int argc, char** argv)
{
    if (argc == 3)
    {
        t_limit = atoi(argv[1]);
        seed = atoi(argv[2]);
    }
    else
    {
        return 0;
    }
    return 1;
}
/*****************************************************************************/
/*****************            Main Function           ************************/
/*****************************************************************************/
int main(int argc, char** argv)
{
    if (GetParam(argc, argv) == 0)
        return 0;
    Initializing();
    if (N >= 500) {
        InitColor();
        for (int i = 0; i < N; ++i)
        {
            cout << Color[i] << endl;
        }
        return 1;
    }
    int p = K;
    while (K)
    {
        clock_t start, finish;
        double totaltime;
        start = clock();

        InitColor();
        InitAdjcaent_Color_Table();
        InitTenure();
        InitConflictNumber();
        int count = 0;
        int breakcount = 0;
        while (f_best)
        {
            Nstep step;
            int delta = NeighborMotion(step, count);
            if (step.v == -1)
                continue;
            Color[step.v] = step.dc;
            TabuTenure[step.v][step.sc] = count + f + (rand() % 10);//*delta

            count++;
            if (f_best <= f)
            {
                ++breakcount;
                if (Max_Non_Improv_One < breakcount) //Termination
                {
                    break;
                }
            }
            else
            {
                breakcount = 0;
            }
            Update(step);   //update Adjcaent_Color_Table & remember the best solution
        }
        if (Max_Non_Improv_One < breakcount)
        {
            if (K == p) {
                for (int i = 0; i < N; ++i)
                {
                    cout << Best_Color[i] << endl;
                }
                return 1;
            }
            else {
                break;
            }
        }
        for (int i = 0; i < N; ++i)
        {
            Global_Best_Color[i] = Best_Color[i];
        }
        finish = clock();
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        --K;
        if (totaltime > t_limit) {
            for (int i = 0; i < N; ++i)
            {
                cout << Global_Best_Color[i] << endl;
            }
            return 1;
        }
    }
    for (int i = 0; i < N; ++i)
    {
        cout << Global_Best_Color[i] << endl;
    }
#ifdef SHOWDEBUG
    cout << "Vertices  Color" << endl;
    for (int i = 0; i < N; ++i)
    {
        cout << i << "---" << Global_Best_Color[i] << "\t";
        if (i % 5 == 4)
        {
            cout << endl;
        }
    }
#endif
    return 1;
}

