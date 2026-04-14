#ifndef LWPART_REFINER_H
#define LWPART_REFINER_H

#include "utils/Types.h"
#include "utils/Configuration.h"
#include "datastructures/Hypergraph.h"
#include "partitioning/Partition.h"
#include "partitioning/PartitionConstraints.h"

namespace consmlp {

/**
 * @brief Refinement statistics
 */
struct RefinementStats {
    Weight initial_cut;
    Weight final_cut;
    Weight improvement;
    uint32_t num_moves;
    uint32_t num_passes;

    RefinementStats()
        : initial_cut(0), final_cut(0), improvement(0)
        , num_moves(0), num_passes(0)
    {}
};

/**
 * @brief Abstract base class for refinement algorithms
 */
class Refiner {
public:
    /**
     * @brief Constructor
     * @param config Configuration parameters
     */
    explicit Refiner(const Configuration& config);

    virtual ~Refiner() = default;

    /**
     * @brief Refine a partition
     * @param hg Hypergraph
     * @param partition Partition to refine (modified in-place)
     * @param constraints Partition constraints
     * @return Refinement statistics
     */
    virtual RefinementStats refine(const Hypergraph& hg,
                                   Partition& partition,
                                   const PartitionConstraints& constraints) = 0;

    virtual const char* getName() const = 0;

protected:
    Configuration config_;
};

/**
 * @brief Greedy FM refinement - optimized version that only initializes cut nodes
 *
 * Key optimizations:
 * - Only initialize gains for nodes connected to cut nets (boundary nodes)
 * - Lazily add adjacent nodes to gain bucket after moves
 * - Adaptive early-quit based on graph size and pass number
 * - Significantly faster than standard FM for large graphs
 */
class GreedyFMRefiner : public Refiner {
public:
    explicit GreedyFMRefiner(const Configuration& config);

    RefinementStats refine(const Hypergraph& hg,
                          Partition& partition,
                          const PartitionConstraints& constraints) override;

    const char* getName() const override { return "GreedyFMRefiner"; }

private:
    /**
     * @brief Gain bucket for O(1) max-gain node retrieval
     */
    class GainBucket {
    public:
        GainBucket();

        /**
         * @brief Initialize with gain range
         * @param min_gain Minimum possible gain
         * @param max_gain Maximum possible gain
         */
        void initialize(int min_gain, int max_gain);

        /**
         * @brief Insert node with gain value
         * @param node_id Node ID
         * @param gain Gain value
         */
        void insert(NodeID node_id, int gain);

        void remove(NodeID node_id);

        /**
         * @brief Update gain for a node already in the bucket
         * @param node_id Node ID
         * @param new_gain New gain value
         */
        void updateGain(NodeID node_id, int new_gain);

        /**
         * @brief Get node with maximum gain (non-destructive peek)
         * @param node_id Output: node ID
         * @param gain Output: gain value
         * @return True if found
         */
        bool getMax(NodeID& node_id, int& gain);

        bool empty() const { return size_ == 0; }
        void clear();

        bool contains(NodeID node_id) const {
            return node_id < node_in_bucket_.size() && node_in_bucket_[node_id];
        }

    private:
        int min_gain_;
        int max_gain_;
        int max_current_gain_;
        std::vector<std::vector<NodeID>> buckets_;
        std::vector<int> node_gains_;
        std::vector<bool> node_in_bucket_;
        size_t size_;

        inline size_t getBucketIndex(int gain) const {
            return static_cast<size_t>(gain - min_gain_);
        }
    };

    /**
     * @brief Move record for rollback to best-cut checkpoint
     */
    struct MoveRecord {
        NodeID node_id;
        PartitionID from_partition;
        PartitionID to_partition;
        Weight cut_delta;

        MoveRecord(NodeID n, PartitionID f, PartitionID t, Weight d)
            : node_id(n), from_partition(f), to_partition(t), cut_delta(d) {}
    };

    /**
     * @brief Perform one greedy FM pass with adaptive early-quit
     * @param hg Hypergraph
     * @param partition Current partition
     * @param constraints Capacity constraints
     * @param pass_number Current pass index (0-indexed)
     * @return Cut improvement achieved in this pass
     */
    Weight performPass(const Hypergraph& hg,
                      Partition& partition,
                      const PartitionConstraints& constraints,
                      uint32_t pass_number);

    /**
     * @brief Initialize gains only for boundary nodes (nodes on cut nets)
     * @param hg Hypergraph
     * @param partition Current partition
     */
    void initializeBoundaryGains(const Hypergraph& hg,
                                 const Partition& partition);

    /**
     * @brief Compute FM gain for moving node to target partition
     * @param hg Hypergraph
     * @param partition Current partition
     * @param node_id Node to evaluate
     * @param to_partition Target partition
     * @return Gain value (positive = cut reduction)
     */
    int computeGain(const Hypergraph& hg,
                    const Partition& partition,
                    NodeID node_id,
                    PartitionID to_partition) const;

    /**
     * @brief Update neighbor gains and lazily add new boundary nodes to bucket
     * @param hg Hypergraph
     * @param partition Current partition (after move)
     * @param moved_node The node that was just moved
     */
    void updateNeighborGainsAndAddToBucket(const Hypergraph& hg,
                                           const Partition& partition,
                                           NodeID moved_node);

    /**
     * @brief Check if moving node to partition is feasible (capacity constraints)
     * @param hg Hypergraph
     * @param partition Current partition
     * @param constraints Capacity constraints
     * @param node_id Node to move
     * @param to_partition Target partition
     * @return True if move is allowed
     */
    bool isValidMove(const Hypergraph& hg,
                    const Partition& partition,
                    const PartitionConstraints& constraints,
                    NodeID node_id,
                    PartitionID to_partition) const;

    std::vector<GainBucket> gain_buckets_;
    std::vector<bool> locked_;
    std::vector<bool> in_bucket_;
};

} // namespace consmlp

#endif // LWPART_REFINER_H
