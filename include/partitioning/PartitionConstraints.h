#ifndef LWPART_PARTITION_CONSTRAINTS_H
#define LWPART_PARTITION_CONSTRAINTS_H

#include "utils/Types.h"
#include "utils/Configuration.h"
#include "datastructures/NodeType.h"
#include "datastructures/Hypergraph.h"
#include "Partition.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <array>

namespace consmlp {

/**
 * @brief Capacity constraints per node type per partition
 */
struct CapacityConstraint {
    Weight min_capacity;
    Weight max_capacity;

    CapacityConstraint() : min_capacity(0), max_capacity(0) {}
    CapacityConstraint(Weight min_cap, Weight max_cap)
        : min_capacity(min_cap), max_capacity(max_cap) {}
};

/**
 * @brief Manages partition constraints (capacity per type per partition)
 */
class PartitionConstraints {
public:
    /**
     * @brief Constructor
     * @param num_partitions Number of partitions
     * @param config Configuration
     */
    PartitionConstraints(PartitionID num_partitions,
                        const Configuration& config);

    /**
     * @brief Initialize balanced constraints (total weight equally divided)
     * @param hg Hypergraph
     * @param imbalance_factor Allowed imbalance (e.g., 0.03 for 3%)
     */
    void initializeBalanced(const Hypergraph& hg, double imbalance_factor);

    /**
     * @brief Initialize per-type balanced constraints
     * DSP/BRAM/IO get relaxed_multiplier * imbalance; IO gets an additional 3x
     * @param hg Hypergraph
     * @param imbalance_factor Base imbalance factor
     * @param relaxed_multiplier Multiplier for DSP/BRAM/IO (default: 3.0)
     */
    void initializeBalancedWithTypes(const Hypergraph& hg,
                                     double imbalance_factor,
                                     double relaxed_multiplier = 3.0);

    /**
     * @brief Initialize constraints for one step of sequential (1 vs k-1) recursive split
     *
     * Used in ConsMLP_lw's non-power-of-2 recursive path.
     * Partition 0 gets (left_count/total_count) of resources,
     * Partition 1 gets (right_count/total_count) of resources.
     *
     * @param left_count  Final partitions routed to partition 0 (left sub-problem)
     * @param right_count Final partitions routed to partition 1 (right sub-problem)
     * @param total_count left_count + right_count
     * @param total_weight_by_type Per-type total weight of the current sub-hypergraph
     * @param imbalance   Base imbalance factor
     * @param use_types   true = apply per-type; false = constrain only LUT (total weight)
     * @param relaxed_multiplier Capacity relaxation for DSP/BRAM/IO (default 3.0)
     */
    void initializeForBipartition(PartitionID left_count,
                                  PartitionID right_count,
                                  PartitionID total_count,
                                  const std::array<Weight, NUM_NODE_TYPES>& total_weight_by_type,
                                  double imbalance,
                                  bool use_types,
                                  double relaxed_multiplier = 3.0);

    /**
     * @brief Initialize constraints from XML constraint file (absolute capacities)
     * @param xml_filename Path to the XML constraint file
     * @return Number of partitions parsed
     */
    PartitionID initializeFromXML(const std::string& xml_filename);

    inline bool isXMLConstraintMode() const { return xml_constraint_mode_; }
    inline void setXmlConstraintMode(bool enabled) { xml_constraint_mode_ = enabled; }

    /**
     * @brief Set capacity constraint for a specific type and partition
     * @param partition_id Partition ID
     * @param type Node type
     * @param min_capacity Minimum capacity
     * @param max_capacity Maximum capacity
     */
    void setCapacity(PartitionID partition_id, NodeType type,
                    Weight min_capacity, Weight max_capacity);

    /**
     * @brief Get capacity constraint
     * @param partition_id Partition ID
     * @param type Node type
     * @return Capacity constraint struct
     */
    CapacityConstraint getCapacity(PartitionID partition_id,
                                   NodeType type) const;

    inline void addFixedNode(NodeID node_id, PartitionID partition_id) {
        fixed_nodes_[node_id] = partition_id;
    }

    inline bool isNodeFixed(NodeID node_id) const {
        return fixed_nodes_.find(node_id) != fixed_nodes_.end();
    }

    inline PartitionID getFixedPartition(NodeID node_id) const {
        auto it = fixed_nodes_.find(node_id);
        return it != fixed_nodes_.end() ? it->second : INVALID_PARTITION;
    }

    /**
     * @brief Check if adding a node (single type) would violate capacity
     * @param partition_id Target partition
     * @param node_type Node type
     * @param node_weight Weight to add
     * @param current_weight Current partition weight for this type
     */
    bool wouldViolateCapacity(PartitionID partition_id, NodeType node_type,
                             Weight node_weight, Weight current_weight) const;

    /**
     * @brief Check if adding a multi-type node would violate any capacity constraint
     * @param partition_id Target partition
     * @param type_weights Per-type weights of the node
     * @param partition Current partition (for current weight lookup)
     */
    bool wouldViolateCapacityMultiType(PartitionID partition_id,
                                       const TypeWeights& type_weights,
                                       const Partition& partition) const;

    bool isBalanced(const Partition& partition, const Hypergraph& hg) const;

    inline PartitionID getNumPartitions() const { return num_partitions_; }

    /**
     * @brief Check whether any partition has a finite (bounded) max capacity for a type.
     * @param type Node type to inspect
     * @return true if at least one partition constrains this type
     */
    bool hasFiniteCapacityForType(NodeType type) const;

    void printConstraintSummary(const Hypergraph& hg) const;
    void printConstraintViolations(const Partition& partition,
                                   const Hypergraph& hg) const;

private:
    PartitionID num_partitions_;
    std::vector<std::vector<CapacityConstraint>> capacity_constraints_;
    std::unordered_map<NodeID, PartitionID> fixed_nodes_;
    bool xml_constraint_mode_ = false;

    static constexpr size_t NUM_TYPES = NUM_NODE_TYPES;
};

} // namespace consmlp

#endif // LWPART_PARTITION_CONSTRAINTS_H
