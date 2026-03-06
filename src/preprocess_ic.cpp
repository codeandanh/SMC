// preprocess_ic.cpp
// IC preprocess for mygraph.cpp (tinyGraph adjacency-list)
// Chốt yêu cầu:
// 1) No self-loop: skip u==v
// 2) node wht ALWAYS random in (0,1)
// 3) dedup keep-first
// 4) trim lines before processing
//
// Output binary format MUST match tinyGraph::read_bin() in mygraph.cpp:
// [node_id n][double preprocessTime]
// for i=0..n-1:
//   [double node_wht][size_t deg]
//   deg times: [node_id nei][double w]

#include "mygraph.cpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <cctype>

using namespace std;
using namespace mygraph;

// Pack (u,v) into uint64 (u,v <= 2^32-1)
static inline std::uint64_t pack_uv(std::uint32_t u, std::uint32_t v) {
    return (static_cast<std::uint64_t>(u) << 32) | static_cast<std::uint64_t>(v);
}
static inline std::uint32_t unpack_u(std::uint64_t key) {
    return static_cast<std::uint32_t>(key >> 32);
}
static inline std::uint32_t unpack_v(std::uint64_t key) {
    return static_cast<std::uint32_t>(key & 0xFFFFFFFFULL);
}

static inline double clamp_nonneg(double x) { return (x >= 0.0 ? x : 0.0); }
static inline double clamp01(double p) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return 1.0;
    return p;
}

// trim (both ends)
static inline std::string trim_copy(const std::string& s) {
    size_t i = 0, j = s.size();
    while (i < j && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    while (j > i && std::isspace(static_cast<unsigned char>(s[j - 1]))) --j;
    return s.substr(i, j - i);
}

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog
        << " edges.txt output.bin [input_undirected(0/1)=1] [seed=42]\n\n"
        << "edges.txt format:\n"
        << "  - Moi dong: u v [w]\n"
        << "    + u, v: so nguyen khong am (ID goc), co the khong lien tiep\n"
        << "    + w   : optional, mac dinh 1.0\n"
        << "    + w am: clamp ve 0\n"
        << "  - Dong rong / comment (#) se bi bo qua (sau khi TRIM).\n\n"
        << "Rules:\n"
        << "  - Renumber theo THU TU XUAT HIEN LAN DAU (deterministic theo file order)\n"
        << "  - No self-loop: bo qua canh u==v\n"
        << "  - Dedup (u,v): keep-first\n"
        << "  - Normalize incoming: sum_in[v]>0 => w(u->v)/=sum_in[v], clamp [0,1]\n"
        << "  - Node wht ALWAYS random in (0,1)\n"
        << "  - Output stored as directed outgoing adjacency (if input_undirected=1, add both dirs)\n";
}

// Write binary in the format expected by tinyGraph::read_bin()
static bool write_bin_for_tinygraph_readbin(const tinyGraph& g, const std::string& out) {
    std::ofstream fout(out.c_str(), std::ios::out | std::ios::binary);
    if (!fout) return false;

    node_id n = static_cast<node_id>(g.n);
    double preprocessTime = g.preprocessTime;

    fout.write(reinterpret_cast<const char*>(&n), sizeof(node_id));
    fout.write(reinterpret_cast<const char*>(&preprocessTime), sizeof(double));

    for (unsigned i = 0; i < g.n; ++i) {
        const double node_wht = g.adjList[i].wht;
        const size_t deg = g.adjList[i].neis.size();

        fout.write(reinterpret_cast<const char*>(&node_wht), sizeof(double));
        fout.write(reinterpret_cast<const char*>(&deg), sizeof(size_t));

        for (size_t j = 0; j < deg; ++j) {
            const node_id nei_id = g.adjList[i].neis[j].getId();
            const double w = g.adjList[i].neis[j].weight;
            fout.write(reinterpret_cast<const char*>(&nei_id), sizeof(node_id));
            fout.write(reinterpret_cast<const char*>(&w), sizeof(double));
        }
    }
    fout.close();
    return true;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string edges_txt  = argv[1];
    const std::string output_bin = argv[2];

    bool input_undirected = true;
    unsigned seed = 42;

    if (argc >= 4) input_undirected = (std::atoi(argv[3]) != 0);
    if (argc >= 5) seed = static_cast<unsigned>(std::strtoul(argv[4], nullptr, 10));

    std::cout << "Preprocess IC for mygraph.cpp (ALWAYS normalize incoming):\n";
    std::cout << "  edges.txt        = " << edges_txt  << "\n";
    std::cout << "  output.bin       = " << output_bin << "\n";
    std::cout << "  input_undirected = " << (input_undirected ? "true" : "false") << "\n";
    std::cout << "  no self-loop     = true\n";
    std::cout << "  node wht random  = true (0,1)\n";
    std::cout << "  dedup            = keep-first\n";
    std::cout << "  trim             = true\n";
    std::cout << "  seed             = " << seed << "\n";

    std::ifstream fin(edges_txt);
    if (!fin) {
        std::cerr << "Error: cannot open edges file: " << edges_txt << "\n";
        return 1;
    }

    // Renumber by first appearance: orig_id (uint64) -> node_id (uint32)
    std::unordered_map<std::uint64_t, node_id> idmap;
    idmap.reserve(1 << 20);

    node_id next_id = 0;
    auto get_new_id = [&](std::uint64_t orig) -> node_id {
        auto it = idmap.find(orig);
        if (it != idmap.end()) return it->second;
        node_id nid = next_id++;
        idmap.emplace(orig, nid);
        return nid;
    };

    // keep-first dedup for directed edges
    std::unordered_map<std::uint64_t, double> keep_first;
    keep_first.reserve(1 << 20);

    size_t num_lines = 0;
    size_t num_kept_directed = 0;
    size_t num_skipped_directed = 0;
    size_t num_skipped_selfloop = 0;

    auto try_insert_keep_first = [&](node_id a, node_id b, double w_raw) {
        const std::uint64_t key = pack_uv(static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(b));
        if (keep_first.find(key) != keep_first.end()) {
            ++num_skipped_directed;
            return;
        }
        keep_first.emplace(key, w_raw);
        ++num_kept_directed;
    };

    std::string line;
    while (std::getline(fin, line)) {
        std::string t = trim_copy(line);
        if (t.empty() || t[0] == '#') continue;

        ++num_lines;

        std::stringstream ss(t);
        std::uint64_t u0, v0;
        if (!(ss >> u0 >> v0)) continue;

        node_id u = get_new_id(u0);
        node_id v = get_new_id(v0);

        double w = 1.0;
        if (ss >> w) w = clamp_nonneg(w);
        else w = 1.0;

        // No self-loop
        if (u == v) {
            ++num_skipped_selfloop;
            continue;
        }

        // Output stored as directed outgoing adjacency
        try_insert_keep_first(u, v, w);
        if (input_undirected) {
            if (v != u) try_insert_keep_first(v, u, w); // redundant guard
        }
    }
    fin.close();

    const unsigned n = static_cast<unsigned>(next_id);
    if (n == 0) {
        std::cerr << "Error: empty graph (no nodes found).\n";
        return 1;
    }

    std::cout << "  Renumbered nodes n   = " << n << "\n";
    std::cout << "  lines_read           = " << num_lines << "\n";
    std::cout << "  kept_directed        = " << num_kept_directed << "\n";
    std::cout << "  skipped_duplicates   = " << num_skipped_directed << "\n";
    std::cout << "  skipped_selfloops    = " << num_skipped_selfloop << "\n";

    // Deterministic materialization order: sort keys by (u,v)
    std::vector<std::uint64_t> keys;
    keys.reserve(keep_first.size());
    for (const auto& kv : keep_first) keys.push_back(kv.first);

    std::sort(keys.begin(), keys.end(), [](std::uint64_t a, std::uint64_t b) {
        const std::uint32_t au = unpack_u(a), av = unpack_v(a);
        const std::uint32_t bu = unpack_u(b), bv = unpack_v(b);
        if (au != bu) return au < bu;
        return av < bv;
    });

    // Build tinyGraph as directed outgoing adjacency
    tinyGraph g;
    g.n = n;
    g.m = 0;                // sum out-deg (half-edges)
    g.preprocessTime = 0.0;
    g.init_empty_graph();

    // Node weights ALWAYS random in (0,1)
    // Use open interval via eps margins.
    {
        std::mt19937 rng(seed);
        const double eps = 1e-9;
        std::uniform_real_distribution<double> dist(eps, 1.0 - eps);
        for (unsigned i = 0; i < g.n; ++i) g.adjList[i].wht = dist(rng);
    }

    // Compute incoming sums on dedup edge set
    std::vector<double> sum_in(g.n, 0.0);
    for (std::uint64_t key : keys) {
        const node_id v = static_cast<node_id>(unpack_v(key));
        double w_raw = keep_first[key];
        if (w_raw > 0.0 && v < g.n) sum_in[static_cast<size_t>(v)] += w_raw;
    }

    // Insert edges with normalized weights (and clamp [0,1])
    for (std::uint64_t key : keys) {
        const node_id u = static_cast<node_id>(unpack_u(key));
        const node_id v = static_cast<node_id>(unpack_v(key));

        double w_raw = keep_first[key];
        if (w_raw < 0.0) w_raw = 0.0;

        double w = w_raw;
        if (v < g.n) {
            const double s = sum_in[static_cast<size_t>(v)];
            if (s > 0.0) w = w_raw / s;
            else w = 0.0;
        } else {
            w = 0.0;
        }
        w = clamp01(w);

        // directed: only add u -> v
        // add_edge_half also rejects self-loop, but self-loop already filtered.
        std::vector<tinyEdge>::iterator edgeAdded;
        g.add_edge_half(u, v, edgeAdded, w);
    }

    // Ensure sorted lists (safe)
    for (unsigned i = 0; i < g.n; ++i) {
        std::sort(g.adjList[i].neis.begin(), g.adjList[i].neis.end(), tinyEdgeCompare);
    }

    // Write binary compatible with tinyGraph::read_bin()
    if (!write_bin_for_tinygraph_readbin(g, output_bin)) {
        std::cerr << "Error: cannot write binary file: " << output_bin << "\n";
        return 1;
    }

    std::cout << "Done. Saved binary graph to: " << output_bin << "\n";
    std::cout << "Summary:\n";
    std::cout << "  n = " << g.n << "\n";
    std::cout << "  m (sum out-deg) = " << g.m << "\n";
    std::cout << "  normalize_incoming = true\n";

    // Verify by reading back
    std::cout << "Verify read_bin...\n";
    tinyGraph h;
    h.read_bin(output_bin);
    std::cout << "  read_bin: n = " << h.n << ", m(sum deg) = " << h.m << "\n";

    return 0;
}
