#ifndef LWPART_COARSENER_H
#define LWPART_COARSENER_H

#include "utils/Types.h"
#include "utils/Configuration.h"
#include "utils/Profiler.h"
#include "datastructures/Hypergraph.h"
#include "datastructures/HypergraphHierarchy.h"
#include <vector>
#include <limits>

namespace consmlp {

/**
 * @brief Coarsening statistics
 */
struct CoarseningStats {
    NodeID original_nodes;
    NodeID coarse_nodes;
    EdgeID original_nets;
    EdgeID coarse_nets;
    double contraction_ratio;
    uint32_t num_matched_pairs;
    uint32_t num_singletons;

    CoarseningStats()
        : original_nodes(0), coarse_nodes(0)
        , original_nets(0), coarse_nets(0)
        , contraction_ratio(0.0)
        , num_matched_pairs(0), num_singletons(0)
    {}
};

/**
 * @brief Abstract base class for coarsening algorithms
 */
class Coarsener {
public:
    /**
     * @brief Constructor
     * @param config Configuration parameters
     */
    explicit Coarsener(const Configuration& config);

    virtual ~Coarsener() = default;

    /**
     * @brief Coarsen a hypergraph and add to hierarchy
     * @param hierarchy Hypergraph hierarchy
     * @param level_idx Current level index (fine level)
     * @param profiler Optional profiler
     * @return Coarsening statistics
     */
    virtual CoarseningStats coarsen(HypergraphHierarchy& hierarchy,
                                   uint32_t level_idx,
                                   Profiler* profiler = nullptr) = 0;

    virtual const char* getName() const = 0;

    /**
     * @brief Check if should stop coarsening
     * @param num_nodes Number of nodes in current graph
     * @return True if should stop
     */
    bool shouldStopCoarsening(NodeID num_nodes) const;

protected:
    Configuration config_;
    uint32_t curr_level_idx_;

    /**
     * @brief Check if node can be coarsened (IO and ignored nodes cannot)
     * @param hg Hypergraph
     * @param node_id Node ID
     */
    bool canCoarsenNode(const Hypergraph& hg, NodeID node_id) const;

    /**
     * @brief Check if two nodes can be matched together
     * @param hg Hypergraph
     * @param node1 First node
     * @param node2 Second node
     */
    bool canMatchNodes(const Hypergraph& hg, NodeID node1, NodeID node2) const;
};

/**
 * @brief Cluster Matching coarsening using computeClusteringV2
 *
 * Algorithm:
 * - Allow clustering of 2+ nodes into a single coarse node
 * - V2 accumulates contributions from multiple neighbors in the same cluster
 * - Selects cluster join vs new pair based on accumulated contribution comparison
 * - Stop when all nodes matched or contraction ratio > 3.0
 */
class ClusterMatching : public Coarsener {
public:
    explicit ClusterMatching(const Configuration& config);

    CoarseningStats coarsen(HypergraphHierarchy& hierarchy,
                           uint32_t level_idx,
                           Profiler* profiler = nullptr) override;

    const char* getName() const override { return "ClusterMatching"; }

private:
    /**
     * @brief Compute clustering V2 - accumulates cluster contribution
     *
     * Unlike V1, this version accumulates connection contributions for nodes
     * already in the same cluster. For example:
     * - Node a has neighbors b, c, d; b is unclustered, c and d are in cluster C
     * - a-b contribution = 1.5, a-c = 1, a-d = 1
     * - V1 would match a with b (1.5 > 1); V2 calculates cluster C = 1+1=2, so a joins C
     *
     * @param hg Hypergraph
     * @param cluster_id Output: node -> cluster_id mapping
     * @return Number of clusters formed
     */
    NodeID computeClusteringV2(const Hypergraph& hg,
                               std::vector<NodeID>& cluster_id);

    /**
     * @brief Judge if a node can be added to an existing cluster
     * @param hg Hypergraph
     * @param cluster_size Current size of the cluster
     * @param cluster_weight Current total weight of the cluster
     * @param node_weight Weight of node to add
     * @return True if node can join the cluster
     */
    bool adjMatchAreaJudging(const Hypergraph& hg,
                             size_t cluster_size,
                             Weight cluster_weight,
                             Weight node_weight) const;

    /**
     * @brief Judge if two unclustered nodes can be matched (weight check)
     * @param hg Hypergraph
     * @param node1 First node
     * @param node2 Second node
     * @param max_weight Maximum combined weight allowed
     * @return True if nodes can be matched
     */
    bool adjMatchAreaJudging(const Hypergraph& hg,
                             NodeID node1, NodeID node2,
                             Weight max_weight) const;

    /**
     * @brief Build coarse nets with parallel net merging and duplicate pin removal
     * @param fine_hg Fine hypergraph
     * @param node_mapping Mapping from fine to coarse nodes
     * @param coarse_hg Coarse hypergraph (output)
     */
    void buildCoarseNets(const Hypergraph& fine_hg,
                        const std::vector<NodeID>& node_mapping,
                        Hypergraph& coarse_hg);
};

} // namespace consmlp

#endif // LWPART_COARSENER_H
