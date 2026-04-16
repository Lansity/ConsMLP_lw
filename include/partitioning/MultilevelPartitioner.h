#ifndef LWPART_MULTILEVEL_PARTITIONER_H
#define LWPART_MULTILEVEL_PARTITIONER_H

#include "utils/HgrParser.h"
#include "utils/Timer.h"
#include "utils/Configuration.h"
#include "datastructures/Hypergraph.h"
#include "datastructures/HypergraphHierarchy.h"
#include "coarsening/Coarsener.h"
#include "coarsening/MultilevelCoarsener.h"
#include "partitioning/Partition.h"
#include "partitioning/PartitionConstraints.h"
#include "partitioning/PartitionMetrics.h"
#include "partitioning/Partitioner.h"
#include "refinement/Refiner.h"

#include <string>
#include <memory>
#include <mutex>
#include <vector>

namespace consmlp {

/**
 * @brief Per-bipartition statistics
 */
struct BipartitionStats {
    int depth = 0;
    PartitionID part_start = 0;
    PartitionID part_end = 0;
    NodeID num_nodes = 0;
    EdgeID num_nets = 0;
    int num_levels = 0;
    double coarsen_time = 0.0;
    double initial_time = 0.0;
    double refine_time = 0.0;
    Weight cut_size = 0;
};

/**
 * @brief Application class for lightweight multilevel hypergraph partitioning
 *
 * Supports:
 * 1. Direct bipartition (k=2, mode=direct)
 * 2. Direct k-way partitioning (k>2, mode=direct)
 * 3. Recursive bipartition with sequential split (k arbitrary, mode=recursive)
 *    - Uses 1 vs k-1 split at each level for non-power-of-2 k support
 */
class MultilevelPartitionerApp {
public:
    MultilevelPartitionerApp() = default;
    ~MultilevelPartitionerApp() = default;

    bool parseArguments(int argc, char* argv[]);
    int run();

    const Configuration& getConfig() const { return config_; }
    const std::string& getCoarsenAlgo() const { return coarsen_algo_; }
    const std::string& getRefineAlgo() const { return refine_algo_; }
    const std::string& getInitMode() const { return init_mode_; }
    bool useTypeConstraints() const { return use_type_constraints_; }
    bool useXMLConstraints() const { return use_xml_constraints_; }
    const std::string& getXMLFile() const { return xml_file_; }
    double getRelaxedMultiplier() const { return relaxed_multiplier_; }

private:
    std::string input_file_;
    std::string output_file_;
    std::string type_file_;
    std::string xml_file_;
    std::string coarsen_algo_;
    std::string refine_algo_;
    std::string partition_mode_;
    std::string init_mode_;
    double relaxed_multiplier_ = 3.0;
    Configuration config_;

    bool use_type_constraints_ = false;
    bool use_xml_constraints_ = false;
    bool use_recursive_ = false;

    std::unique_ptr<PartitionConstraints> xml_constraints_;
    PartitionID xml_total_partitions_ = 0;
    std::vector<PartitionID> xml_active_partition_ids_;
    double xml_target_utilization_ = 0.60;

    void setLargeNetThreshold(Hypergraph& hg);
    void adjustXMLActivePartitions(const Hypergraph& hg);
    void printConfiguration() const;
    Hypergraph parseHypergraph();

    int runRecursivePartitioning(Hypergraph& hg);
    int runDirectPartitioning(Hypergraph& hg, double parse_time);

    Partition runInitialPartitioning(const HypergraphHierarchy& hierarchy,
                                    const int coarsest_level,
                                    const PartitionConstraints& coarsest_constraints,
                                    bool use_type_constraints,
                                    double relaxed_multiplier,
                                    const std::string& init_mode,
                                    int& best_trial_level);
};

// Utility functions
bool isPowerOfTwo(PartitionID value);

/**
 * @brief Create coarsener (ConsMLP_lw: always ClusterMatching)
 * @param algo Algorithm name (ignored, only cluster supported)
 * @param config Configuration
 */
std::unique_ptr<Coarsener> createCoarsener(const std::string& algo,
                                           const Configuration& config);

/**
 * @brief Create refiner (ConsMLP_lw: always GreedyFMRefiner)
 * @param algo Algorithm name (ignored, only gfm supported)
 * @param config Configuration
 */
std::unique_ptr<Refiner> createRefiner(const std::string& algo,
                                       const Configuration& config);

/**
 * @brief Build a sub-hypergraph induced by a subset of nodes
 * @param parent Original hypergraph
 * @param subset Node IDs from parent to include
 * @return Sub-hypergraph with remapped node IDs [0, subset.size())
 */
Hypergraph buildSubHypergraph(const Hypergraph& parent,
                              const std::vector<NodeID>& subset);

/**
 * @brief Recursive bipartition using sequential (1 vs k-1) splitting
 *
 * At each call, the subset is split into left (1 final partition) and
 * right (target_parts-1 final partitions) using proportional capacity
 * constraints: left = total_weight*(1/target_parts)*(1±imb),
 *              right = total_weight*((k-1)/target_parts)*(1±imb).
 *
 * @param original_hg Original finest hypergraph
 * @param subset Current node subset (original node IDs)
 * @param base_partition_id Starting partition ID for this sub-problem
 * @param target_parts Number of final partitions to produce
 * @param base_config Configuration
 * @param coarsen_algo Coarsening algorithm
 * @param refine_algo Refinement algorithm
 * @param init_mode Initial partitioning mode
 * @param use_type_constraints Enable per-type resource constraints
 * @param relaxed_multiplier Imbalance relaxation for DSP/BRAM/IO
 * @param final_partition Output partition (written to with mutex)
 * @param all_stats Output statistics vector
 * @param xml_constraints XML absolute capacity constraints (nullptr = proportional)
 * @param partition_mutex Mutex protecting final_partition writes
 * @param stats_mutex Mutex protecting all_stats writes
 * @param depth Current recursion depth
 */
void recursiveBipartitionDirect(const Hypergraph& original_hg,
                                const std::vector<NodeID>& subset,
                                PartitionID base_partition_id,
                                PartitionID target_parts,
                                const Configuration& base_config,
                                const std::string& coarsen_algo,
                                const std::string& refine_algo,
                                const std::string& init_mode,
                                bool use_type_constraints,
                                double relaxed_multiplier,
                                Partition& final_partition,
                                std::vector<BipartitionStats>& all_stats,
                                const PartitionConstraints* xml_constraints,
                                std::mutex& partition_mutex,
                                std::mutex& stats_mutex,
                                int depth);

void printUsage(const char* prog_name);
void printHypergraphInfo(const Hypergraph& hg, bool use_types = false);
void printCoarseningResults(const MultilevelCoarsener& coarsener,
                            const HypergraphHierarchy& hierarchy,
                            double time_s);
void printRefinementResults(const std::vector<RefinementStats>& stats,
                           double time_s);
void printFinalResults(const Hypergraph& hg, const Partition& partition,
                      const Configuration& config);

} // namespace consmlp

#endif // LWPART_MULTILEVEL_PARTITIONER_H
