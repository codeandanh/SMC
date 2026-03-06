#ifndef ALGS_CPP
#define ALGS_CPP

#include "mygraph.cpp"
#include "obj.cpp"
#include <set>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <cfloat>
#include <unordered_map>
#include <malloc.h>
#include <limits>

using namespace std;
using namespace mygraph;

enum Algs
{
   ALG1,
   ALG2,
   ALG3,
   ALG4,// thuật toán cũ trong bài csonet 2023
   GREEDY,
   SINGLE,
   MULTI,
   USM
};

resultsHandler allResults;
mutex mtx;

class Args
{
public:
   Algs alg;
   string graphFileName;
   string outputFileName = "";
   size_t k = 2;
   tinyGraph g;
   double tElapsed;
   double wallTime;
   Logger logg;
   double epsi = 0.1;
   double delta = 0.1;
   double c = 1;
   size_t N = 1;
   size_t P = 10;
   bool plusplus = false;
   double tradeoff = 0.5;
   bool quiet = false;
   bool lazy = false;
   size_t nThreads = 1;
   double buget;
   double T = 1;
   Args() {}
   Args(const Args &args) : alg(args.alg),
                            graphFileName(args.graphFileName),
                            outputFileName(args.outputFileName),
                            k(args.k),
                            buget(args.buget),
                            T(args.T),
                            g(args.g),
                            epsi(args.epsi),
                            delta(args.delta),
                            c(args.c),
                            N(args.N),
                            P(args.P),
                            plusplus(args.plusplus),
                            lazy(args.lazy),
                            nThreads(args.nThreads)
   {
   }
};

class MyPair
{
public:
   node_id u;
   double gain;

   MyPair() {}
   MyPair(node_id a,
          double g)
   {
      u = a;
      gain = g;
   }

   MyPair(const MyPair &rhs)
   {
      u = rhs.u;
      gain = rhs.gain;
   }

   void operator=(const MyPair &rhs)
   {
      u = rhs.u;
      gain = rhs.gain;
   }
};

struct gainLT
{
   bool operator()(const MyPair &p1, const MyPair &p2)
   {
      return p1.gain < p2.gain;
   }
} gainLTobj;

struct revgainLT
{
   bool operator()(const MyPair &p1, const MyPair &p2)
   {
      return (p1.gain > p2.gain);
   }
} revgainLTobj;

void reportResults(size_t nEvals, size_t obj, size_t maxMem = 0)
{
   mtx.lock();
   allResults.add("obj", obj);
   allResults.add("nEvals", nEvals);
   allResults.add("mem", maxMem);
   mtx.unlock();
}

void random_set(tinyGraph &g, vector<bool> &C, vector<size_t> &A)
{
   C.assign(g.n, false);
   double prob = 0.5;
   uniform_real_distribution<double> unidist(0, 1);
   for (size_t i = 0; i < A.size(); ++i)
   {
      if (unidist(gen) < prob)
      {
         C[A[i]] = true;
      }
   }
}

class Alg1
{
   random_device rd;
   mt19937 gen;
   double epsi, T;
   Args &myArgs;
   size_t k, nEvals = 0;
   tinyGraph &g;
   double alpha, beta;
   double c_min = 1.0;
   double c_max = 99.0;
   double delta = 0.1;
   int number_pass = 0;
   VectorHash hasher;

public:
   Alg1(Args &args) : gen(rd()), myArgs(args), g(args.g)
   {
      k = args.k;
      epsi = args.epsi;
      T = args.T;
   }
   void print(Solution sol)
   {
      cout << "u=" << sol.u << " j=" << sol.j << " curren_f=" << sol.curren_f << " s size: " << sol.s.size() << endl;
   }
   vector<Solution> get_U(double lower_bound, double c_min, double upper_bound)
   {
      vector<Solution> U;
      int j = 1;
      double value;
      while (true)
      {
         value = c_min * pow(1 + epsi, j - 1);
         if (value < lower_bound)
         {
            j++;
            continue;
         }
         if (value > upper_bound)
            break;
         Solution sol;
         sol.u = value;
         sol.j = j;
         U.push_back(sol);
         j++;
      }
      return U;
   }
   void get_c_min_max(double &c_min, double &c_max)
   {
      c_max = -1;
      c_min = numeric_limits<double>::max();
//#pragma omp parallel for
      for (node_id i = 0; i < g.n; i++)
      {
         if (g.adjList[i].wht > c_max)
         {
//#pragma omp critical
            {
               c_max = g.adjList[i].wht;
            }
         }
         if (g.adjList[i].wht < c_min)
         {
//#pragma omp critical
            {
               c_min = g.adjList[i].wht;
            }
         }
      }
      //cout<<"c_min: "<<c_min<<" c_max: "<<c_max<<endl;
   }
   double get_cost(vector<node_id> S)
   {
      double c = 0;
      //#pragma omp parallel for reduction(+ : c)
      for (int i = 0; i < S.size(); i++)
      {
         c += g.adjList[S[i]].wht;
      }
      return c;
   }
   double run_alg1()
   {
      get_c_min_max(c_min, c_max);
      nEvals = 0;
      vector<Solution> U = get_U(c_min / (1 + epsi), c_min, c_max * g.n);
      alpha = 2 * (1 - epsi);
      beta = (1 + epsi) * alpha / epsi;
      Solution best_sol;
      double min_cost = std::numeric_limits<double>::max();
      bool flag = false;
      for (node_id e = 0; e < g.n; e++)
      {
         std::unordered_map<int, double> myMap;
         for (Solution &sol : U)
         {
            if (sol.curren_f < (1 - epsi) * T)
            {
               if (sol.cost + g.adjList[e].wht <= beta * sol.u)
               {
                  sol.s.push_back(e);
                  int hashValue = hasher(sol.s);
                  sol.s.pop_back();
                  auto it = myMap.find(hashValue);
                  if (it == myMap.end())
                  {
                     sol.s.push_back(e);
                     double tmp = compute_valSet(nEvals, g, sol.s);
                     myMap[hashValue] = tmp;
                     sol.s.pop_back();
                  }
                  double tmp = myMap[hashValue];
                  if (tmp - sol.curren_f >= g.adjList[e].wht * alpha * T / (beta * sol.u))
                  {
                     sol.cost += g.adjList[e].wht;
                     sol.curren_f = tmp;
                     sol.s.push_back(e);
                     if (sol.curren_f >= (1 - epsi) * T && sol.cost < min_cost)
                     {
                        flag = true;
                        min_cost = sol.cost;
                        best_sol = sol;
                     }
                  }
               }
               if (g.adjList[e].wht <= (1 + epsi) * sol.u)
               {
                  vector<node_id> tmp;
                  tmp.push_back(e);
                  double tmp_e = compute_valSet(nEvals, g, tmp);
                  if (tmp_e > sol.f_e)
                  {
                     sol.e = e;
                     sol.f_e = tmp_e;
                     sol.cost_e = g.adjList[e].wht;
                     if (sol.f_e >= (1 - epsi) * T && sol.cost_e < min_cost)
                     {
                        min_cost = sol.cost;
                        best_sol = sol;
                        flag = false;
                     }
                  }
               }
            }
         }
      }
      number_pass = 1;
      cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      if (flag)
         cout << "ALG1," << T << "," << best_sol.curren_f << "," << nEvals << "," << best_sol.cost << "," << number_pass;
      else
         cout << "ALG1," << T << "," << best_sol.f_e << "," << nEvals << "," << best_sol.cost_e << "," << number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "," << totalMemoryUsed << "," << totalMemoryUnused;
      return best_sol.curren_f;
   }
   void get_C_min_max_prime(double &c_min_prime, double &c_max_prime)
   {
      std::multiset<std::pair<node_id, double>, CustomCompare> mySet;
      int count = 0;
      for (node_id e = 0; e < g.n; e++)
      {
         mySet.insert({e, g.adjList[e].wht});
         count++;
      }
      vector<node_id> v;
      node_id uj;
      for (auto it = mySet.begin(); it != mySet.end(); it++)
      {
         v.push_back(it->first);
         double tmp = compute_valSet(nEvals, g, v);
         if (tmp >= T)
         {
            uj = it->first;
            c_min_prime = epsi * it->second / g.n;
            c_max_prime = (it->first + 1) * it->second;
            break;
         }
      }
   }
   vector<node_id> get_V0(double const c_min_prime, double const c_max_prime, double &C, double &c_V0)
   {
      c_V0 = 0;
      C = 0;
      vector<node_id> v;
      for (node_id e = 0; e < g.n; e++)
      {
         if (g.adjList[e].wht < c_min_prime)
         {
            v.push_back(e);
            c_V0 += g.adjList[e].wht;
         }
         else if (g.adjList[e].wht <= c_max_prime && g.adjList[e].wht >= c_min_prime)
         {
            C += g.adjList[e].wht;
         }
      }
      return v;
   }
   double run_alg2()
   {
      get_c_min_max(c_min, c_max);
      nEvals = 0;
      alpha = 2 * (1 - epsi);
      beta = (1 + epsi) * alpha / (2 - alpha);
      double c_min_prime, c_max_prime;
      // first pass
      get_C_min_max_prime(c_min_prime, c_max_prime);
      number_pass++;
      vector<node_id> V0;
      double C;
      double cost_V0;
      // second pass
      V0 = get_V0(c_min_prime, c_max_prime, C, cost_V0);
      number_pass++;
      double f_V0 = compute_valSet(nEvals, g, V0);
      // double T_prime=T-f_V0;
      vector<Solution> U = get_U(c_min_prime / (1 + epsi), c_min_prime, C);
      //#pragma omp parallel for
      for (Solution &sol : U)
      {
         sol.s = V0;
         sol.cost = cost_V0;
         sol.curren_f = f_V0;
      }
      // third pass
      Solution best_solution;
      double min_cost = numeric_limits<double>::max();
      bool flag;
      for (node_id e = 0; e < g.n; e++)
      {
         if (g.adjList[e].wht >= c_min_prime && g.adjList[e].wht <= c_max_prime)
         {
            std::unordered_map<int, double> myMap;
            for (Solution &sol : U)
            {
               if (sol.curren_f < (1 - epsi) * T)
               {
                  if (sol.cost + g.adjList[e].wht <= beta * sol.u && sol.cost + g.adjList[e].wht<min_cost)
                  {
                     sol.s.push_back(e);
                     int hashValue = hasher(sol.s);
                     sol.s.pop_back();
                     auto it = myMap.find(hashValue);
                     if (it == myMap.end())
                     {
                        sol.s.push_back(e);
                        double tmp = compute_valSet(nEvals, g, sol.s);
                        myMap[hashValue] = tmp;
                        sol.s.pop_back();
                     }
                     double tmp_f = myMap[hashValue];

                     if (tmp_f - sol.curren_f >= g.adjList[e].wht * alpha * T / (beta * sol.u))
                     {
                        sol.s.push_back(e);
                        sol.curren_f = tmp_f;
                        sol.cost += g.adjList[e].wht;
                        if (sol.curren_f >= (1 - epsi) * T && sol.cost < min_cost)
                        {

                           min_cost = sol.cost;
                           best_solution = sol;
                           flag = true;
                        }
                     }
                  }
                  if (sol.cost_e <= (1 + epsi) * sol.u)
                  {
                     vector<node_id> tmp = V0;
                     tmp.push_back(e);
                     double tmp_fe = compute_valSet(nEvals, g, tmp);
                     if (tmp_fe > sol.f_e)
                     {
                        sol.f_e = tmp_fe;
                        sol.cost_e = g.adjList[e].wht;
                        if (sol.f_e >= (1 - epsi) * T)
                        {
                           min_cost = sol.cost_e;
                           best_solution = sol;
                           flag = false;
                        }
                     }
                  }
               }
            }
         }
      }
      number_pass++;
      cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      if (flag)
         cout << "ALG2," << T << "," << best_solution.curren_f << "," << nEvals << "," << best_solution.cost << "," << number_pass;
      else
         cout << "ALG2_e," << T << "," << best_solution.f_e << "," << nEvals << "," << best_solution.cost_e << "," << number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "," << totalMemoryUsed << "," << totalMemoryUnused;
      return 0.0;
   }
   double run_alg3()
   {
      number_pass = 0;
      get_c_min_max(c_min, c_max);
      nEvals = 0;
      alpha = 2 * (1 - epsi);
      beta = (1 + epsi) * alpha / (2 - alpha);
      double c_min_prime, c_max_prime;
      // first pass
      get_C_min_max_prime(c_min_prime, c_max_prime);
      number_pass++;
      vector<node_id> V0;
      double C;
      double cost_V0;
      // second pass
      V0 = get_V0(c_min_prime, c_max_prime, C, cost_V0);
      number_pass++;
      // cout<<"V0 size: "<<V0.size()<<" C: "<<C<<" cost_V0: "<<cost_V0<<endl;
      double f_V0 = compute_valSet(nEvals, g, V0);

      vector<Solution> U = get_U(c_min_prime / (1 + epsi), c_min_prime, C);
      //#pragma omp parallel for
      for (Solution &sol : U)
      {
         sol.s = V0;
         sol.cost = cost_V0;
         sol.curren_f = f_V0;
         sol.check = vector<bool>(g.n, false);
      }
      // Multi-pass
      Solution best_solution;
      double min_cost = numeric_limits<double>::max();
      bool flag;
      int ell = ceil((1.0 / epsi) * log(1.0 / delta));
      int counter = 0;
      bool check = false;
      // cout<<"ell: "<<ell<<endl;
      for (int i = 0; i < ell; i++)
      {
         number_pass++;
         check = false;
         for (node_id e = 0; e < g.n; e++)
         {
            if (g.adjList[e].wht >= c_min_prime && g.adjList[e].wht <= c_max_prime)
            {
               std::unordered_map<int, double> myMap;
               for (Solution &sol : U)
               {
                  if (sol.check[e] == false && sol.curren_f < (1 - delta) * T)
                  {
                     if (sol.cost + g.adjList[e].wht <= beta * sol.u && sol.cost + g.adjList[e].wht <=min_cost)
                     {
                        sol.s.push_back(e);
                        int hashValue = hasher(sol.s);
                        sol.s.pop_back();
                        auto it = myMap.find(hashValue);
                        if (it == myMap.end())
                        {
                           sol.s.push_back(e);
                           double tmp = compute_valSet(nEvals, g, sol.s);
                           myMap[hashValue] = tmp;
                           sol.s.pop_back();
                        }
                        double tmp_f = myMap[hashValue];

                        if (tmp_f - sol.curren_f >= g.adjList[e].wht * (1 - epsi) * (T - sol.curren_f) / sol.u)
                        {
                           sol.s.push_back(e);
                           sol.check[e] = true;
                           sol.curren_f = tmp_f;
                           sol.cost += g.adjList[e].wht;
                           if (sol.curren_f >= (1 - epsi) * T && sol.cost < min_cost)
                           {
                              min_cost = sol.cost;
                              best_solution = sol;
                              check = true;
                           }
                        }
                     }
                  }
               }
            }
         }
         // cout << "------ ALG3," << T << "," << best_solution.curren_f << "," << nEvals << "," << best_solution.cost<<endl;
         if (check == false)
            counter++;
         if (counter == 2)
            break;
      }
      cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      cout << "ALG3," << T << "," << best_solution.curren_f << "," << nEvals << "," << best_solution.cost << "," << number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "," << totalMemoryUsed << "," << totalMemoryUnused;
      return 0.0;
   }
   int Findj(vector<myNode> A)
   {
      int j=0;
      vector<node_id> VJ;
      double f_Vj;
      while (true)
      {
         if (j >= g.n) break;
         VJ.push_back(A[j].u);
         f_Vj = compute_valSet(nEvals, g, VJ);
         if (f_Vj >= T) break;
         j++;
      }
      return j;
   }
   void Divide_V(double c_min, double c_max, vector<node_id> &V0,double &c_V0,vector<node_id> &VV,double &c_VV, double &c_0)
   {
      c_0 = std::numeric_limits<double>::max();
      c_V0=0;
      for (node_id u = 0; u < g.n; u++)
      {
         if (g.adjList[u].wht < c_min)
         {
            V0.push_back(u);
            c_V0+=g.adjList[u].wht;
            continue;
         }
         else if (g.adjList[u].wht <= c_max)
         {
            VV.push_back(u);
            c_0 = min(c_0, g.adjList[u].wht);
            c_VV += g.adjList[u].wht;
         }
      }
   }
   vector<int> Get_U(double c_0,double c_VV)
   {
      int i = 0;
      vector<int> U;
      while (true)
      {
         double tmp = pow((1 + epsi), i);
         if (tmp < c_0)
         {
            i++;
            continue;
         }
         if (tmp > c_VV)
            break;
         U.push_back(i);
         i++;
      }
      return U;
   }
   double run_alg4()// thuật toán cũ trong bài csonet 2023
   {
      nEvals = 0;
      alpha = 2 * (1 - epsi);
      beta = (1 + epsi) * alpha / (2 - alpha);
      vector<myNode> A;
      for (node_id u = 0; u < g.n; u++)
      {
         A.push_back(myNode(u, g.adjList[u].wht));
      }
      std::sort(A.begin(), A.end(), compareByCost);
      //cout<<"A size: "<<A.size()<<endl;
      //tìm j
      int j = Findj(A);
      double c_min = epsi * g.adjList[A[j].u].wht / g.n;
      double c_max = (j + 1) * g.adjList[A[j].u].wht;
      //cout<<"j: "<<j<<" c_min: "<<c_min<<" c_max: "<<c_max<<endl;
      vector<node_id> V0, VV;
      double c_0 ,c_VV,c_v0;
      Divide_V(c_min,c_max,V0,c_v0,VV,c_VV,c_0);
      //cout<<"V0 size: "<<V0.size()<<" VV size: "<<VV.size()<<endl;
      double f_V0 = compute_valSet(nEvals, g, V0);
      vector<int> U=Get_U(c_0,c_VV);
      //cout<<"U size: "<<U.size()<<endl;
      vector<node_id> S,SS;
      double f_current=0;
      double c_s=0;
      double f_max=0,cost_min=std::numeric_limits<double>::max();
      int number_pass=1;
      for (int v = 0; v < U.size(); v++)
      {
         S=V0;
         c_s=c_v0;
         node_id emax;
         double f_emax=f_V0,cost_eMax;
         f_current=f_V0;
         for (int i = 0; i < VV.size(); i++)
         {
            node_id e=VV[i];
            if (c_s + g.adjList[e].wht <= beta * pow((1 + epsi), v))
            {
               S.push_back(e);
               double f_tmp = compute_valSet(nEvals, g, S);
               S.pop_back();

               if ((f_tmp - f_current) / g.adjList[e].wht >= alpha * T / (beta * pow((1 + epsi), v)))
               {
                  S.push_back(e);
                  c_s += g.adjList[e].wht;
                  f_current = f_tmp;
               }
            }
            if (g.adjList[e].wht <= pow((1 + epsi), v+1))
            {
               V0.push_back(e);
               double f_e = compute_valSet(nEvals, g, V0);
               V0.pop_back();
               f_emax = max(f_emax,f_e);
               cost_eMax=g.adjList[e].wht;
            }
            if(f_current>=alpha * T /2 || f_emax>=alpha * T /2) break;
         }
         number_pass++;
         if(f_current>=alpha * T /2 && c_s<cost_min)
         {
            cost_min=c_s;
            f_max=f_current;
         }

         if(f_emax >= alpha * T /2 && cost_eMax<cost_min)
         {
            cost_min=cost_eMax;
            f_max=f_emax;
         }        
         if(f_max >= alpha * T /2) break;
      }
      //cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      cout << "ALG4," << T << "," << f_max << "," << nEvals << "," << cost_min<<","<<number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout<<","<<totalMemoryUsed<<","<<totalMemoryUnused;
      return f_current;
   }
};
class Multi
{
   random_device rd;
   mt19937 gen;
   vector<node_id> S, A, B;
   double epsi, T;
   int i = 1;
   Args &myArgs;
   size_t k, nEvals = 0;
   tinyGraph &g;
   int number_pass = 0;
   Logger logg;

public:
   Multi(Args &args) : gen(rd()), myArgs(args), g(args.g)
   {
      k = args.k;
      epsi = args.epsi;
      T = args.T;
   }
   double get_cost(vector<node_id> S)
   {
      double c = 0;
      for (int i = 0; i < S.size(); i++)
      {
         c += g.adjList[S[i]].wht;
      }
      return c;
   }
   double USM(vector<node_id> S, double &cost)
   {
      double f=compute_valSet(nEvals, g, S);
      cost = get_cost(S);
      return f;
   }
   double stream(double opt, double &cost)
   {
      int l = 2.0 / epsi;
      vector<vector<node_id>> S(l + 1, vector<node_id>());
      vector<double> f_S(l + 1, 0.0);
      vector<double> c_S(l + 1, 0.0);
      vector<bool> check_S0(g.n, false);
      for (int u = 0; u < g.n; u++)
      {
         if (g.adjList[u].wht > opt) continue;
         int jj=-1;
         double tmp_f=0;
         for (int j = 1; j <= l; j++)
         {
            S[j].push_back(u);
            tmp_f = compute_valSet(nEvals, g, S[j]);
            S[j].pop_back();

            if ((tmp_f - f_S[j]) / g.adjList[u].wht >= epsi * T / (2 * opt))
            {
               jj=j;
               break;
               
            }
         } 
         if(jj!=-1)
         {
            S[jj].push_back(u);
            f_S[jj]=tmp_f;
            c_S[jj] += g.adjList[u].wht;
               
            S[0].push_back(u);
            c_S[0]+= g.adjList[u].wht;
            if (c_S[jj] > 2 * opt / epsi) break;
         }      
      }
      number_pass++;
      double f_max=0, cost_max=0;
      for(int j = 1; j <= l; j++)
      {
         if(f_S[j]>f_max)
         {
            f_max=f_S[j];
            cost_max=c_S[j];
         }
      }
      double cost_usm=0;
      double fmax_usm = USM(S[0], cost_usm);
      if(f_max > fmax_usm)
      {
         cost = cost_max;
         return f_max;
      }
      else
      {
         cost = cost_usm;
         return fmax_usm;
      }
      return 0.0;
   }
   double get_min_cost()
   {
      double opt = numeric_limits<double>::max();
      for (int i = 0; i < g.n; i++)
      {
         if (g.adjList[i].wht < opt)
            opt = g.adjList[i].wht;
      }
      return opt;
   }
   double get_total_cost()
   {
      double total_cost = 0.0;
      for (int i = 0; i < g.n; i++)
      {
         total_cost += g.adjList[i].wht;
      }
      return total_cost;
   }
   double run_multi()
   {
      number_pass=0;
      nEvals = 0;
      double gamma = 1.0, f = 0;
      double opt = get_min_cost();
      //cout<<"min cost: "<<opt<<endl;
      double size = 0;
      while (true)
      {
         f = stream(opt, size);
         if (f >= gamma * (1 - epsi) * T)
            break;
         opt = (1 + epsi) * opt;
         //cout<<"opt: "<<opt<<" f:"<<f<<endl;
      }
      cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      cout << "Multi," << T << "," << f << "," << nEvals << "," << size << "," << number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "," << totalMemoryUsed << "," << totalMemoryUnused;
      return 0.0;
   }
   double compute_val_u(node_id u)
   {
      vector<node_id> s_tmp;
      s_tmp.push_back(u);
      return compute_valSet(nEvals, g, s_tmp);
   }
   double run_single()
   {
      nEvals = 0;
      number_pass=0;
      double b = get_total_cost(), gamma = 1.0;
      int h = (2.0 / epsi), n = g.n;
      double L = -1;
      vector<int> sigma;
      {
         int j = 0;
         while (true)
         {
            double tmp = pow((1 + epsi), j);
            if (tmp > b)
               break;
            sigma.push_back(j);
            j++;
         }
      }
      int size_sigma=sigma.size();
      vector<node_id> s_tmp;
      vector<vector<node_id>> htmp(h + 1, s_tmp);
      vector<vector<vector<node_id>>> S(size_sigma, htmp);
      vector<vector<double>> f_S(size_sigma, vector<double>(h + 1, 0));
      vector<vector<double>> cost_s(size_sigma, vector<double>(h + 1, 0));
      vector<double> max_f_sig(size_sigma, 0);
      vector<double> size_sig(size_sigma, 0);

      double fmax = 0;
      int size_max = 0;
      int index_sig;

      for (int u = 0; u < g.n; u++)
      {
         double tmp_fu = compute_val_u(u);
         if (tmp_fu / g.adjList[u].wht > epsi * T / (2 * L))
            L = epsi * T * g.adjList[u].wht / (2 * tmp_fu);
         for (int t : sigma)
         {
            double sigg = pow((1 + epsi), t);
            if (sigg < L || sigg > b) continue;

            bool check=true;
            for (int i = 1; i <= h; i++)
            {
               if(cost_s[t][i] >= 2 * sigg / epsi)
               {
                  check=false;
                  break;
               }
            }
            int ii = -1;
            double tmp_f = 0;
            for (int i = 1; i <= h; i++)
            {
               if (check==true && g.adjList[u].wht <= sigg)
               {
                  S[t][i].push_back(u);
                  tmp_f = compute_valSet(nEvals, g, S[t][i]);
                  S[t][i].pop_back();
                  if (tmp_f - f_S[t][i] >= g.adjList[u].wht * epsi * T / (2 * sigg))
                  {
                     ii=i;
                     break;                   
                  }
               }
            }
            if(ii!=-1)
            {
               S[t][ii].push_back(u);
               f_S[t][ii] = tmp_f;
               cost_s[t][ii] += g.adjList[u].wht;

               S[t][0].push_back(u);
               double size;
               f_S[t][0] = USM(S[t][0], size);
               cost_s[t][0] = size;
            }
            double max_f=-1;
            int max_i=-1;
            for(int i = 0; i <= h; i++ )
            {
               if(f_S[t][i]>max_f)
               {
                  max_f=f_S[t][i];
                  max_i=i;
               }
            }
            if (max_f >= gamma * (1 - epsi) * T)
            {
               b = sigg;
               index_sig = t;
            }
         }
      }
      number_pass++;
      double max_f = 0;
      double cost_min = numeric_limits<double>::max();
      for (int i = 0; i <= h; i++)
      {
         if (f_S[index_sig][i] >=(1-epsi)*T && cost_s[index_sig][i] < cost_min)
         {
            max_f = f_S[index_sig][i];
            cost_min = cost_s[index_sig][i];
         }
      }
      cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      cout << "single," << T << "," << max_f << "," << nEvals << "," << cost_min << "," << number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "," << totalMemoryUsed << "," << totalMemoryUnused;
      return 0.0;
   }
   double runUSM()
   {
      number_pass = 0;
      nEvals = 0;
      vector<node_id> S;
      for (node_id u = 0; u < g.n; u++)
      {
         S.push_back(u);
      }
      double c;
      double f = USM(S, c);
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "alg,f,q,cost,pass,memory,memoryUnused,time" << endl;
      cout << std::setprecision(15) << "f:" << f << ", cost:" << c << ", Quries: " << nEvals << ", Memory:" << totalMemoryUsed << ",MemoryUnused:" << totalMemoryUnused;
      return f;
   }
};

class Greedy
{
   random_device rd;
   mt19937 gen;
   double epsi, T;
   Args &myArgs;
   size_t k, nEvals = 0;
   tinyGraph &g;
   int number_pass = 0;

public:
   Greedy(Args &args) : gen(rd()), myArgs(args), g(args.g)
   {
      k = args.k;
      epsi = args.epsi;
      T = args.T;
   }
   double run()
   {
      number_pass = 0;
      nEvals = 0;
      vector<node_id> S;
      vector<bool> check_s(g.n, false);
      double f_s = 0;
      double cost_s = 0;
      while (true)
      {
         double max_f = 0, delta = 0;
         int u_max = -1;
         for (node_id u = 0; u <= g.n; u++)
         {
            if (check_s[u])
               continue;

            S.push_back(u);
            double tmp_f = compute_valSet(nEvals, g, S);
            S.pop_back();

            double tmp_delta = (min(tmp_f, T) - f_s) / g.adjList[u].wht;
            if (tmp_delta > delta)
            {
               delta = tmp_delta;
               u_max = u;
               max_f = tmp_f;
            }
         }
         number_pass++;
         if (u_max == -1)
            break;
         S.push_back(u_max);
         check_s[u_max] = true;
         f_s = max_f;
         cost_s += g.adjList[u_max].wht;
         if (f_s >= T * (1 - epsi))
            break;
      }
      cout << "alg,T,f,q,cost,pass,memory,memoryUnused,time" << endl;
      cout << "Greedy," << T << "," << f_s << "," << nEvals << "," << cost_s << "," << number_pass;
      struct mallinfo info = mallinfo();
      int totalMemoryUsed = info.uordblks;
      int totalMemoryUnused = info.fordblks;
      cout << "," << totalMemoryUsed << "," << totalMemoryUnused;
      return f_s;
   }
};
#endif
