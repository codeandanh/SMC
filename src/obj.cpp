#ifndef OBJ_CPP
#define OBJ_CPP

#include <omp.h>
#include "mygraph.cpp"

// (IC needs these; harmless for others)
#include <vector>
#include <utility>
#include <random>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <cmath>

using namespace std;
using namespace mygraph;

vector<bool> emptySetVector;

#ifdef MAXCOV
bool monotone = true;
size_t marge(size_t &nEvals, tinyGraph &g, node_id u, vector<bool> &set,
             vector<bool> &cov)
{

   if (set[u])
      return 0;

   ++nEvals;

   return g.getDegreeMinusSet(u, cov) + 1;
}

size_t marge(size_t &nEvals, tinyGraph &g, node_id u, vector<bool> &set)
{
   if (set[u])
      return 0;

   vector<bool> cov(g.n, false);
   size_t val = 0;
   for (node_id u = 0; u < g.n; ++u)
   {
      if (set[u])
      {
         if (!cov[u])
         {
            cov[u] = true;
            val += 1;
         }
         vector<tinyEdge> &neis = g.adjList[u].neis;
         for (size_t j = 0; j < neis.size(); ++j)
         {
            node_id v = neis[j].target;
            if (!cov[v])
            {
               cov[v] = true;
               val += 1;
            }
         }
      }
   }

   return marge(nEvals, g, u, set, cov);
}

size_t compute_valSet(size_t &nEvals, tinyGraph &g, vector<bool> &set,
                      vector<bool> &cov = emptySetVector)
{
   ++nEvals;
   cov.assign(g.n, false);
   size_t val = 0;
   for (node_id u = 0; u < g.n; ++u)
   {
      if (set[u])
      {
         if (!cov[u])
         {
            cov[u] = true;
            val += 1;
         }
         vector<tinyEdge> &neis = g.adjList[u].neis;
         for (size_t j = 0; j < neis.size(); ++j)
         {
            node_id v = neis[j].target;
            if (!cov[v])
            {
               cov[v] = true;
               val += 1;
            }
         }
      }
   }

   return val;
}

size_t compute_valSet(size_t &nEvals, tinyGraph &g, vector<node_id> &set)
{
   ++nEvals;
   vector<bool> cov(g.n, false);

   size_t val = 0;
   for (size_t i = 0; i < set.size(); ++i)
   {
      node_id u = set[i];
      if (!cov[u])
      {
         cov[u] = true;
         val += 1;
      }
      vector<tinyEdge> &neis = g.adjList[u].neis;
      for (size_t j = 0; j < neis.size(); ++j)
      {
         node_id v = neis[j].target;
         if (!cov[v])
         {
            cov[v] = true;
            val += 1;
         }
      }
   }

   return val;
}
#endif

#ifdef REVMAX_MON
bool monotone = true;
vector<double> alpha;

void init_alpha(tinyGraph &g)
{
   uniform_real_distribution<double> unidist(0, 1);
   alpha.assign(g.n, 0.0);
   mt19937 gen(0); // same sequence each time

   for (node_id u = 0; u < g.n; ++u)
   {
      alpha[u] = unidist(gen);
   }
}

double compute_valSet(size_t &nEvals, tinyGraph &g, vector<bool> &set,
                      vector<bool> &cov = emptySetVector)
{
   if (alpha.size() == 0)
   {
      init_alpha(g);
   }

   ++nEvals;
   cov.assign(g.n, false);
   double val = 0;

   for (node_id u = 0; u < g.n; ++u)
   {
      vector<tinyEdge> &neis = g.adjList[u].neis;
      double valU = 0.0;
      for (size_t j = 0; j < neis.size(); ++j)
      {
         node_id v = neis[j].target;
         if (set[v])
         {
            valU += neis[j].weight;
         }
      }
      valU = pow(valU, alpha[u]);
      val += valU;
   }

   return val;
}

double compute_valSet(size_t &nEvals, tinyGraph &g, vector<node_id> &sset)
{
   if (alpha.size() == 0)
   {
      init_alpha(g);
   }

   vector<bool> set(g.n, false);
   for (size_t i = 0; i < sset.size(); ++i)
   {
      set[sset[i]] = true;
   }

   ++nEvals;

   double val = 0;

   for (node_id u = 0; u < g.n; ++u)
   {
      vector<tinyEdge> &neis = g.adjList[u].neis;
      double valU = 0.0;
      for (size_t j = 0; j < neis.size(); ++j)
      {
         node_id v = neis[j].target;
         if (set[v])
         {
            valU += neis[j].weight;
         }
      }
      valU = pow(valU, alpha[u]);
      val += valU;
   }

   return val;
}

double marge(size_t &nEvals, tinyGraph &g, node_id x, vector<bool> &set,
             vector<bool> &cov = emptySetVector)
{
   if (alpha.size() == 0)
   {
      init_alpha(g);
   }

   if (set[x])
      return 0;

   vector<tinyEdge> &neis = g.adjList[x].neis;
   double gain = 0.0;
   for (size_t j = 0; j < neis.size(); ++j)
   {
      node_id v = neis[j].target;
      vector<tinyEdge> &neisV = g.adjList[v].neis;
      double valV = 0.0;
      double valVwithX = 0.0;
      for (size_t k = 0; k < neisV.size(); ++k)
      {
         node_id w = neisV[k].target;
         if (w != x)
         {
            if (set[w])
            {
               valV += neisV[k].weight;
               valVwithX += neisV[k].weight;
            }
         }
         else
         {
            valVwithX += neisV[k].weight;
         }
      }

      if (valV == 0)
         gain += pow(valVwithX, alpha[v]);
      else
         gain += pow(valVwithX, alpha[v]) - pow(valV, alpha[v]);
   }
   ++nEvals;
   return gain;
}

#endif

#ifdef MAXCUT
bool monotone = false;
signed long marge(size_t &nEvals, tinyGraph &g, node_id u, vector<bool> &set, vector<bool> &ancillary = emptySetVector)
{

   if (set[u])
      return 0;

   ++nEvals;

   signed long m;
   double mx = 2 * g.getWeightedDegreeMinusSet(u, set);
   double my = g.getWeightedDegree(u);

   m = (mx - my);

   return m;
}

size_t compute_valSet(size_t &nEvals, tinyGraph &g, vector<bool> &set)
{
   ++nEvals;
   size_t val = 0;
#pragma omp parallel for reduction(+:val)
   for (node_id u = 0; u < g.n; ++u)
   {
      vector<tinyEdge> &neis = g.adjList[u].neis;
      for (size_t j = 0; j < neis.size(); ++j)
      {
         node_id v = neis[j].target;
         if ((set[u] && !set[v]) || (!set[u] && set[v]))
            val += neis[j].weight;
      }
   }

   return val / 2;
}

double compute_valSet(size_t &nEvals, tinyGraph &g, vector<node_id> &set_id)
{
   vector<bool> set(g.n, false);
   for (size_t i = 0; i < set_id.size(); ++i)
   {
      set[set_id[i]] = true;
   }

   ++nEvals;
   double val = 0;
#pragma omp parallel for reduction(+:val)
   for (node_id u = 0; u < g.n; ++u)
   {
      vector<tinyEdge> &neis = g.adjList[u].neis;
      for (size_t j = 0; j < neis.size(); ++j)
      {
         node_id v = neis[j].target;
         if ((set[u] && !set[v]) || (!set[u] && set[v]))
            val += neis[j].weight;
      }
   }

   return val / 2;
}

#endif

#ifdef REVMAX_NM
bool monotone = false;
vector<double> alpha;

void init_alpha(tinyGraph &g)
{
   uniform_real_distribution<double> unidist(0, 1);
   alpha.assign(g.n, 0.0);
   mt19937 gen(0); // same sequence each time

   for (node_id u = 0; u < g.n; ++u)
   {
      alpha[u] = unidist(gen);
   }
}

double compute_valSet(size_t &nEvals, tinyGraph &g, vector<bool> &set, vector<bool> &cov = emptySetVector)
{
   if (alpha.size() == 0)
   {
      init_alpha(g);
   }

   ++nEvals;
   cov.assign(g.n, false);
   double val = 0;

   for (node_id u = 0; u < g.n; ++u)
   {
      if (!set[u])
      {
         vector<tinyEdge> &neis = g.adjList[u].neis;
         double valU = 0.0;
         for (size_t j = 0; j < neis.size(); ++j)
         {
            node_id v = neis[j].target;
            if (set[v])
            {
               valU += neis[j].weight;
            }
         }
         valU = pow(valU, alpha[u]);
         val += valU;
      }
   }

   return val;
}

double compute_valSet(size_t &nEvals, tinyGraph &g, vector<node_id> &sset)
{
   if (alpha.size() == 0)
   {
      init_alpha(g);
   }
   vector<bool> set(g.n, false);
   for (size_t i = 0; i < sset.size(); ++i)
   {
      set[sset[i]] = true;
   }

   ++nEvals;

   double val = 0;
#pragma omp parallel for reduction(+:val)
   for (node_id u = 0; u < g.n; ++u)
   {
      if (!set[u])
      {
         vector<tinyEdge> &neis = g.adjList[u].neis;
         double valU = 0.0;
         for (size_t j = 0; j < neis.size(); ++j)
         {
            node_id v = neis[j].target;
            if (set[v])
            {
               valU += neis[j].weight;
            }
         }
         valU = sqrt(valU);
         val += valU;
      }
   }

   return val;
}

double marge(size_t &nEvals, tinyGraph &g, node_id x, vector<bool> &set,
             vector<bool> &cov = emptySetVector)
{
   if (alpha.size() == 0)
   {
      init_alpha(g);
   }
   if (set[x])
      return 0;

   double loss = 0.0;

   vector<tinyEdge> &neis = g.adjList[x].neis;
   double valX = 0.0;
   for (size_t j = 0; j < neis.size(); ++j)
   {
      node_id v = neis[j].target;
      if (set[v])
      {
         valX += neis[j].weight;
      }
   }

   valX = pow(valX, alpha[x]);

   loss = valX;

   double gain = 0.0;
   for (size_t j = 0; j < neis.size(); ++j)
   {
      node_id v = neis[j].target;
      vector<tinyEdge> &neisV = g.adjList[v].neis;
      double valV = 0.0;
      double valVwithX = 0.0;
      for (size_t k = 0; k < neisV.size(); ++k)
      {
         node_id w = neisV[k].target;
         if (w != x)
         {
            if (set[w])
            {
               valV += neisV[k].weight;
               valVwithX += neisV[k].weight;
            }
         }
         else
         {
            valVwithX += neisV[k].weight;
         }
      }

      gain += pow(valVwithX, alpha[v]) - pow(valV, alpha[v]);
   }

   ++nEvals;
   return gain - loss;
}

#endif

#ifdef IMG
bool monotone = false;
signed long marge(size_t &nEvals, tinyGraph &g, node_id u, vector<bool> &set, vector<bool> &ancillary = emptySetVector)
{

   if (set[u])
      return 0;

   ++nEvals;

   signed long m;
   double mx = 2 * g.getWeightedDegreeMinusSet(u, set);
   double my = g.getWeightedDegree(u);

   m = (mx - my);

   return m;
}

size_t compute_valSet(size_t &nEvals, tinyGraph &g, vector<bool> &set)
{
   ++nEvals;
   size_t val = 0;
#pragma omp parallel for reduction(+:val)
   for (node_id u = 0; u < g.n; ++u)
   {
      vector<tinyEdge> &neis = g.adjList[u].neis;
      for (size_t j = 0; j < neis.size(); ++j)
      {
         node_id v = neis[j].target;
         if ((set[u] && !set[v]) || (!set[u] && set[v]))
            val += neis[j].weight;
      }
   }

   return val / 2;
}

double compute_valSet(size_t &nEvals, tinyGraph &g, vector<node_id> &set_id)
{
   vector<bool> set(g.n, false);
   for (size_t i = 0; i < set_id.size(); ++i)
   {
      set[set_id[i]] = true;
   }

   ++nEvals;
   double val = 0;

#pragma omp parallel for reduction(+:val)
   for (node_id u = 0; u < g.n; ++u)
   {
      double max = 0;
      for (node_id v = 0; v < g.n; v++)
      {
         if (set[u] && set[v])
         {
            val -= (double)g.getEdgeWeight(u, v) / (double)g.n;
         }
         if (set[u] == false && set[v] == true)
         {
            double tmp = g.getEdgeWeight(u, v);
            if (tmp > max)
               max = tmp;
         }
      }
      val += max;
   }

   return val;
}

#endif

// ======================================================
// IC objective (monotone): f(S) = E[ |A(S)| ]
// env:
//   KIC_MC    (default 100)
//   KIC_SEED  (if set -> reproducible; if NOT set -> truly-random per call)
// Notes per your requirements:
//   - interface must match obj.cpp
//   - marge must be truly random (before/after independent, unless KIC_SEED is set)
//   - marge increments nEvals only once
// ======================================================
#ifdef IC

#ifndef _OPENMP
#error "IC requires OpenMP. Please compile with -fopenmp (and link with OpenMP)."
#endif

bool monotone = true;

namespace detail_ic_obj {

static inline double clamp01(double p) {
   if (p <= 0.0) return 0.0;
   if (p >= 1.0) return 1.0;
   return p;
}

static inline bool env_exists(const char* name) {
   return (std::getenv(name) != nullptr);
}

static inline std::size_t read_env_size_t(const char* name, std::size_t defv) {
   if (const char* s = std::getenv(name)) {
      char* end = nullptr;
      unsigned long long v = std::strtoull(s, &end, 10);
      if (end && *end == '\0') return static_cast<std::size_t>(v);
   }
   return defv;
}

static inline std::uint64_t read_env_u64(const char* name, std::uint64_t defv) {
   if (const char* s = std::getenv(name)) {
      char* end = nullptr;
      unsigned long long v = std::strtoull(s, &end, 10);
      if (end && *end == '\0') return static_cast<std::uint64_t>(v);
   }
   return defv;
}

// splitmix64 để seed theo iteration (ổn định theo it)
static inline std::uint64_t splitmix64(std::uint64_t x) {
   x += 0x9e3779b97f4a7c15ULL;
   x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
   x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
   return x ^ (x >> 31);
}

// “ngẫu nhiên theo lần gọi”: random_device + thời gian + địa chỉ
static inline std::uint64_t runtime_seed64() {
   std::random_device rd;
   std::uint64_t r = (static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd());

   const std::uint64_t t =
       static_cast<std::uint64_t>(
           std::chrono::high_resolution_clock::now().time_since_epoch().count());

   std::uint64_t addr_mix = reinterpret_cast<std::uintptr_t>(&rd);
   return splitmix64(r ^ t ^ addr_mix);
}

struct ICParams {
   std::size_t mc;
   std::uint64_t base_seed;
};

static inline ICParams read_params_per_call() {
   ICParams p;
   p.mc = read_env_size_t("KIC_MC", 100);

   // Nếu set KIC_SEED -> reproducible; nếu không -> thật ngẫu nhiên theo mỗi lần gọi
   p.base_seed = env_exists("KIC_SEED")
       ? read_env_u64("KIC_SEED", 42ULL)
       : runtime_seed64();
   return p;
}

// Core evaluation: f(S) = E[|A(S)|]
// IMPORTANT: does NOT touch nEvals; caller controls counting.
static inline double evaluate_core(const mygraph::tinyGraph& g,
                                  const std::vector<bool>& set,
                                  const ICParams& prm)
{
   const std::size_t n = static_cast<std::size_t>(g.n);
   if (n == 0) return 0.0;

   // seeds
   std::vector<mygraph::node_id> seeds;
   seeds.reserve(n);
   for (mygraph::node_id u = 0; u < (mygraph::node_id)g.n; ++u) {
      if (set[(std::size_t)u]) seeds.push_back(u);
   }

   if (prm.mc == 0) return 0.0;

   double sum_spread = 0.0;

#pragma omp parallel
   {
      std::vector<std::uint32_t> seen(n, 0);
      std::uint32_t stamp = 1;

      auto bump_stamp = [](std::uint32_t& st, std::vector<std::uint32_t>& arr) {
         ++st;
         if (st == 0) { // overflow wrap
            std::fill(arr.begin(), arr.end(), 0);
            st = 1;
         }
      };

      std::uniform_real_distribution<double> uni(0.0, 1.0);

#pragma omp for schedule(static) reduction(+:sum_spread)
      for (std::size_t it = 0; it < prm.mc; ++it) {
         // RNG theo it: ổn định theo iteration và không phụ thuộc schedule/num_threads
         std::mt19937_64 rng(splitmix64(prm.base_seed ^ static_cast<std::uint64_t>(it)));

         bump_stamp(stamp, seen);
         std::size_t activated_cnt = 0;

         std::vector<mygraph::node_id> frontier;
         frontier.reserve(seeds.size());

         // init from seeds
         for (mygraph::node_id s : seeds) {
            if ((std::size_t)s >= n) continue;
            if (seen[(std::size_t)s] == stamp) continue;
            seen[(std::size_t)s] = stamp;
            frontier.push_back(s);
            ++activated_cnt;
         }

         // IC BFS-like using adjList[u].neis as outgoing edges
         while (!frontier.empty()) {
            std::vector<mygraph::node_id> nxt;
            nxt.reserve(frontier.size());

            for (mygraph::node_id uu : frontier) {
               const auto& neis = g.adjList[(std::size_t)uu].neis;
               for (std::size_t j = 0; j < neis.size(); ++j) {
                  mygraph::node_id v = neis[j].getId();
                  if ((std::size_t)v >= n) continue;
                  if (seen[(std::size_t)v] == stamp) continue;

                  const double p = clamp01(neis[j].weight);
                  if (p <= 0.0) continue;

                  if (uni(rng) < p) {
                     seen[(std::size_t)v] = stamp;
                     nxt.push_back(v);
                     ++activated_cnt;
                  }
               }
            }
            frontier.swap(nxt);
         }

         sum_spread += static_cast<double>(activated_cnt);
      }
   } // omp parallel

   return sum_spread / static_cast<double>(prm.mc);
}

} // namespace detail_ic_obj

// compute_valSet for vector<bool>
double compute_valSet(size_t &nEvals, tinyGraph &g, vector<bool> &set,
                      vector<bool> &cov /* unused */ = emptySetVector)
{
   (void)cov;
   ++nEvals;

   const auto prm = detail_ic_obj::read_params_per_call();
   return detail_ic_obj::evaluate_core(g, set, prm);
}

// compute_valSet for vector<node_id>
double compute_valSet(size_t &nEvals, tinyGraph &g, vector<node_id> &set_id)
{
   vector<bool> set(g.n, false);
   for (size_t i = 0; i < set_id.size(); ++i)
   {
      node_id u = set_id[i];
      if (u < g.n) set[(size_t)u] = true;
   }

   ++nEvals;

   const auto prm = detail_ic_obj::read_params_per_call();
   return detail_ic_obj::evaluate_core(g, set, prm);
}

// marge = f(S ∪ {x}) - f(S)
// per your requirement: truly random -> before/after independent seeds (unless KIC_SEED is set)
double marge(size_t &nEvals, tinyGraph &g, node_id x, vector<bool> &set,
             vector<bool> &ancillary /* unused */ = emptySetVector)
{
   (void)ancillary;
   if (x >= g.n) return 0.0;
   if (set[(size_t)x]) return 0.0;

   // before: independent seed (if KIC_SEED not set)
   const auto prm1 = detail_ic_obj::read_params_per_call();
   const double before = detail_ic_obj::evaluate_core(g, set, prm1);

   vector<bool> set_after = set;
   set_after[(size_t)x] = true;

   // after: independent seed (if KIC_SEED not set)
   const auto prm2 = detail_ic_obj::read_params_per_call();
   const double after = detail_ic_obj::evaluate_core(g, set_after, prm2);

   ++nEvals; // per your requirement: only +1
   return after - before;
}

#endif // IC

#endif
