#include "coarsening/Coarsener.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <iostream>

namespace consmlp {

namespace {
constexpr Index kMaxRatingNetDegree = 50;
constexpr Index kMaxMatchingNetDegree = 50;
constexpr size_t kMaxCoarseNetSize = 100;  // Reduced for faster coarsening
}

// ========== Coarsener Base Class ==========

Coarsener::Coarsener(const Configuration& config)
    : config_(config)
{
}

bool Coarsener::shouldStopCoarsening(NodeID num_nodes) const {
    return num_nodes <= config_.coarsening_threshold;
}

bool Coarsener::canCoarsenNode(const Hypergraph& hg, NodeID node_id) const {
    // Cannot coarsen ignored nodes
    if (hg.isNodeIgnored(node_id)) {
        return false;
    }
    
    // IO nodes should not be coarsened
    if (hg.getNodeType(node_id) == NodeType::IO) {
        return false;
    }
    
    return true;
}

bool Coarsener::canMatchNodes(const Hypergraph& hg, 
                              NodeID node1, NodeID node2) const {
    // Both nodes must be coarsenable
    if (!canCoarsenNode(hg, node1) || !canCoarsenNode(hg, node2)) {
        return false;
    }
    
    // Fixed nodes: can only match if both fixed to same partition
    bool node1_fixed = hg.isNodeFixed(node1);
    bool node2_fixed = hg.isNodeFixed(node2);
    
    if (node1_fixed || node2_fixed) {
        // If one is fixed and the other is not, cannot match
        if (node1_fixed != node2_fixed) {
            return false;
        }
        // Both are fixed - would need partition info to check
        // For now, allow matching fixed nodes (refinement will handle it)
    }
    
    // NOTE: Different types CAN be matched now
    // Per-type weights are tracked in each coarsened node
    // This allows better coarsening quality while preserving type constraints
    
    return true;
}

ClusterMatching::ClusterMatching(const Configuration& config)
    : Coarsener(config)
{
}

CoarseningStats ClusterMatching::coarsen(HypergraphHierarchy& hierarchy,
                                         uint32_t level_idx) {
    CoarseningStats stats;
    curr_level_idx_ = level_idx;
    HypergraphLevel& fine_level = hierarchy.getLevel(curr_level_idx_);
    const Hypergraph& fine_hg = fine_level.getHypergraph();
    
    stats.original_nodes = fine_hg.getNumNodes();
    stats.original_nets = fine_hg.getNumNets();
    
    // Compute clustering (multi-node matching)
    // V2 considers accumulated cluster contribution for better quality
    std::vector<NodeID> cluster_id(stats.original_nodes, INVALID_NODE);
    NodeID num_clusters = computeClusteringV2(fine_hg, cluster_id);
    
    // Build coarse graph
    HypergraphLevel& coarse_level = hierarchy.addCoarserLevel();
    Hypergraph& coarse_hg = coarse_level.getHypergraph();
    
    // Compute cluster per-type weights (accumulate from all fine nodes in cluster)
    std::vector<TypeWeights> cluster_type_weights(num_clusters);
    for (NodeID i = 0; i < num_clusters; ++i) {
        for (size_t t = 0; t < NUM_NODE_TYPES; ++t) {
            cluster_type_weights[i][t] = 0;
        }
    }
    
    for (NodeID fine_node = 0; fine_node < stats.original_nodes; ++fine_node) {
        NodeID cid = cluster_id[fine_node];
        if (cid != INVALID_NODE) {
            // Accumulate per-type weights from fine node
            const TypeWeights& fine_weights = fine_hg.getNodeTypeWeights(fine_node);
            for (size_t t = 0; t < NUM_NODE_TYPES; ++t) {
                cluster_type_weights[cid][t] += fine_weights[t];
            }
        }
    }
    
    // Add coarse nodes with per-type weights
    for (NodeID i = 0; i < num_clusters; ++i) {
        coarse_hg.addNodeWithTypeWeights(cluster_type_weights[i]);
    }
    
    // Build coarse nets with parallel net merging
    buildCoarseNets(fine_hg, cluster_id, coarse_hg);
    
    coarse_hg.finalize();
    
    stats.coarse_nodes = num_clusters;
    stats.coarse_nets = coarse_hg.getNumNets();
    stats.contraction_ratio = static_cast<double>(stats.original_nodes) / stats.coarse_nodes;
    
    // Count matched pairs (approximate: nodes in clusters of size >= 2)
    std::vector<size_t> cluster_sizes(num_clusters, 0);
    for (NodeID fine_node = 0; fine_node < stats.original_nodes; ++fine_node) {
        NodeID cid = cluster_id[fine_node];
        if (cid != INVALID_NODE) {
            cluster_sizes[cid]++;
        }
    }
    
    stats.num_matched_pairs = 0;
    for (NodeID cid = 0; cid < num_clusters; ++cid) {
        if (cluster_sizes[cid] >= 2) {
            stats.num_matched_pairs += cluster_sizes[cid] - 1;
        }
    }
    stats.num_singletons = stats.coarse_nodes - stats.num_matched_pairs;
    
    // Set up mappings
    for (NodeID fine_node = 0; fine_node < stats.original_nodes; ++fine_node) {
        NodeID coarse_node = cluster_id[fine_node];
        fine_level.setCoarserNode(fine_node, coarse_node);
        fine_level.addFinerNode(coarse_node, fine_node);
    }
    fine_level.finalizeMapping();
    
    return stats;
}

NodeID ClusterMatching::computeClusteringV2(const Hypergraph& hg,
                                            std::vector<NodeID>& cluster_id) {
    constexpr Index kMaxNetDegree = 50;  // Skip large nets for performance
    
    const NodeID num_nodes = hg.getNumNodes();
    
    // Track cluster information
    std::vector<bool> is_clustered(num_nodes, false);
    std::vector<NodeID> cluster_leader(num_nodes, INVALID_NODE);  // Leader of each cluster
    std::vector<size_t> cluster_size(num_nodes, 0);  // Size of cluster led by this node
    std::vector<Weight> cluster_weight(num_nodes, 0);  // Weight of cluster led by this node
    uint32_t max_cluster_size = 12; 
    if (curr_level_idx_ < 3)
        max_cluster_size = 10;
    else if (curr_level_idx_ < 7)
        max_cluster_size = 5;
    else   
        max_cluster_size = 4;
    NodeID num_clusters = 0;
    
    // Visit nodes in ascending order of weight/degree
    std::vector<NodeID> order(num_nodes);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](NodeID a, NodeID b) {
        Weight wa = hg.getNodeWeight(a);
        Weight wb = hg.getNodeWeight(b);
        if (wa == wb) {
            return hg.getNodeDegree(a) < hg.getNodeDegree(b);
        }
        return wa < wb;
    });

    auto total_weight = config_.total_node_weight;
    auto admissible_weight = num_nodes < 300 ? total_weight / 40 : total_weight / 80;

    auto quan_weight_factor = [&](Weight weight) -> double {
        return static_cast<double>(weight) / (static_cast<double>(total_weight) / static_cast<double>(num_nodes));
    };

    double net_scale_factor = 1.0;
    for (NodeID node : order) {
        // Skip if already clustered
        if (is_clustered[node]) continue;
        
        // Skip if cannot coarsen this node
        if (!canCoarsenNode(hg, node)) {
            // Create singleton cluster
            cluster_id[node] = num_clusters;
            cluster_leader[node] = node;
            cluster_size[node] = 1;
            cluster_weight[node] = hg.getNodeWeight(node);
            is_clustered[node] = true;
            num_clusters++;
            continue;
        }
        
        // Key difference from V1: separate tracking for unclustered nodes and clusters
        // unclustered_candidates: neighbor node -> contribution
        // cluster_candidates: cluster leader -> accumulated contribution
        std::unordered_map<NodeID, double> unclustered_candidates;
        std::unordered_map<NodeID, double> cluster_candidates;
        auto nets = hg.getNodeNets(node);
        for (const EdgeID* net_it = nets.first; net_it != nets.second; ++net_it) {
            EdgeID net_id = *net_it;
            Index net_size = hg.getNetSize(net_id);
            // Skip very small or very large nets
            if (net_size <= 1 || net_size > kMaxNetDegree) continue;
            if (num_nodes < 10000) 
                net_scale_factor = 1.0;
            else
                net_scale_factor = static_cast<double>(net_size - 1);
            double contribution = static_cast<double>(hg.getNetWeight(net_id)) / net_scale_factor;
            auto nodes = hg.getNetNodes(net_id);
            for (const NodeID* node_it = nodes.first; node_it != nodes.second; ++node_it) {
                NodeID neighbor = *node_it;
                if (neighbor == node) continue;
                if (!canMatchNodes(hg, node, neighbor)) continue;
                
                
                if (is_clustered[neighbor]) {
                    // Neighbor is in a cluster - accumulate contribution to cluster leader
                    NodeID leader = cluster_leader[neighbor];
                    // double quan_weight = contribution / quan_weight_factor(hg.getNodeWeight(neighbor));
                    double quan_weight = contribution / quan_weight_factor(cluster_weight[leader]);

                    if (leader == INVALID_NODE || leader >= num_nodes) continue;
                    
                    // Check weight constraint
                    if (cluster_weight[leader] + hg.getNodeWeight(node) > admissible_weight) {
                        continue;
                    }
                    
                    // Accumulate contribution to this cluster
                    if (cluster_candidates.find(leader) == cluster_candidates.end()) {
                        cluster_candidates[leader] = quan_weight;
                    }
                    else {
                        cluster_candidates[leader] += quan_weight;
                    }
                } else {
                    // Neighbor is unclustered
                    double quan_weight = contribution / quan_weight_factor(hg.getNodeWeight(neighbor));

                    if (!adjMatchAreaJudging(hg, node, neighbor, admissible_weight)) {
                        continue;
                    }
                    if (unclustered_candidates.find(neighbor) == unclustered_candidates.end()) {
                        unclustered_candidates[neighbor] = quan_weight;
                    }
                    else {
                        unclustered_candidates[neighbor] += quan_weight;
                    }
                }
            }
        }
        
        // Find best unclustered candidate
        NodeID best_unclustered = INVALID_NODE;
        double best_unclustered_score = 0.0;
        for (const auto& p : unclustered_candidates) {
            if (p.second > best_unclustered_score) {
                best_unclustered_score = p.second;
                best_unclustered = p.first;
            }
        }
        
        // Find best cluster candidate
        NodeID best_cluster_leader = INVALID_NODE;
        double best_cluster_score = 0.0;
        for (const auto& p : cluster_candidates) {
            NodeID leader = p.first;
            // Check cluster size constraint
            if (cluster_size[leader] >= max_cluster_size) continue;
            // NOTE: Type matching removed - per-type weights are tracked separately
            
            if (p.second > best_cluster_score) {
                best_cluster_score = p.second;
                best_cluster_leader = leader;
            }
        }
        
        // No candidates found - create singleton
        if (best_unclustered == INVALID_NODE && best_cluster_leader == INVALID_NODE) {
            cluster_id[node] = num_clusters;
            cluster_leader[node] = node;
            cluster_size[node] = 1;
            cluster_weight[node] = hg.getNodeWeight(node);
            is_clustered[node] = true;
            num_clusters++;
            continue;
        }
        
        Weight node_weight = hg.getNodeWeight(node);
        
        // Compare: cluster accumulated score vs unclustered single node score
        if (best_cluster_score >= best_unclustered_score && best_cluster_leader != INVALID_NODE) {
            // Merge into existing cluster (cluster has higher or equal accumulated contribution)
            NodeID leader = best_cluster_leader;
            cluster_id[node] = cluster_id[leader];  // Use leader's cluster ID
            cluster_leader[node] = leader;
            cluster_size[leader]++;
            cluster_weight[leader] += node_weight;
            is_clustered[node] = true;
        } else if (best_unclustered != INVALID_NODE) {
            // Create new cluster with unclustered neighbor
            cluster_id[node] = num_clusters;
            cluster_id[best_unclustered] = num_clusters;
            cluster_leader[node] = node;  // node is the leader
            cluster_leader[best_unclustered] = node;
            cluster_size[node] = 2;
            cluster_weight[node] = node_weight + hg.getNodeWeight(best_unclustered);
            is_clustered[node] = true;
            is_clustered[best_unclustered] = true;
            num_clusters++;
        } else {
            // Fallback: create singleton (should not reach here)
            cluster_id[node] = num_clusters;
            cluster_leader[node] = node;
            cluster_size[node] = 1;
            cluster_weight[node] = node_weight;
            is_clustered[node] = true;
            num_clusters++;
        }
    }
    
    // Renumber clusters to be contiguous
    std::vector<NodeID> cluster_remap(num_clusters, INVALID_NODE);
    NodeID new_cluster_count = 0;
    
    for (NodeID node = 0; node < num_nodes; ++node) {
        NodeID old_cid = cluster_id[node];
        if (old_cid != INVALID_NODE) {
            if (cluster_remap[old_cid] == INVALID_NODE) {
                cluster_remap[old_cid] = new_cluster_count++;
            }
            cluster_id[node] = cluster_remap[old_cid];
        }
    }
    
    return new_cluster_count;
}

/*
 * Reserved for future cluster-area extension.
 * Currently not used by computeClusteringV2.
 *
 * bool ClusterMatching::adjMatchAreaJudging(const Hypergraph& hg,
 *                                           size_t cluster_size,
 *                                           Weight cluster_weight,
 *                                           Weight node_weight) const {
 *     // If cluster has only 2 nodes, always allow joining
 *     if (cluster_size <= 2) {
 *         return true;
 *     }
 *
 *     // For larger clusters, apply constraints
 *     // Limit cluster size to prevent too aggressive clustering
 *     constexpr size_t kMaxClusterSize = 4;
 *     if (cluster_size >= kMaxClusterSize) {
 *         return false;
 *     } else {
 *         return true;
 *     }
 *
 *     // Limit cluster weight (prevent very heavy clusters)
 *     // Use average node weight as reference
 *     Weight avg_weight = hg.getTotalNodeWeight() / hg.getNumNodes();
 *     Weight max_cluster_weight = avg_weight * kMaxClusterSize;
 *
 *     if (cluster_weight + node_weight > max_cluster_weight) {
 *         return false;
 *     }
 *
 *     return true;
 * }
 */

void ClusterMatching::buildCoarseNets(const Hypergraph& fine_hg,
                                      const std::vector<NodeID>& node_mapping,
                                      Hypergraph& coarse_hg) {
    const uint64_t max_coarse_net_size = std::max(config_.large_net_threshold, static_cast<uint64_t>(300));

    // Hash function for vector of sorted node IDs
    struct VectorHash {
        size_t operator()(const std::vector<NodeID>& vec) const {
            size_t hash = vec.size();
            for (NodeID node : vec) {
                hash ^= std::hash<NodeID>()(node) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
            return hash;
        }
    };
    const VectorHash vector_hash;
    
    // Map hash -> candidate net indices (collision-resolved by exact vector compare).
    // This avoids storing a second full vector copy as unordered_map key.
    std::unordered_map<size_t, std::vector<size_t>> hash_to_indices;
    hash_to_indices.reserve(fine_hg.getNumNets() / 2);
    
    std::vector<std::vector<NodeID>> net_keys;
    std::vector<Weight> net_weights;
    net_keys.reserve(fine_hg.getNumNets() / 2);
    net_weights.reserve(fine_hg.getNumNets() / 2);
    
    // Reusable buffer to avoid repeated allocations
    std::vector<NodeID> coarse_nodes_buf;
    coarse_nodes_buf.reserve(max_coarse_net_size);
    
    for (EdgeID fine_net = 0; fine_net < fine_hg.getNumNets(); ++fine_net) {
        coarse_nodes_buf.clear();
        
        auto fine_nodes = fine_hg.getNetNodes(fine_net);
        for (const NodeID* it = fine_nodes.first; it != fine_nodes.second; ++it) {
            NodeID coarse_node = node_mapping[*it];
            if (coarse_node != INVALID_NODE) {
                coarse_nodes_buf.push_back(coarse_node);
            }
        }
        
        if (coarse_nodes_buf.size() <= 1) continue;
        if (coarse_nodes_buf.size() > max_coarse_net_size) continue;

        // Fast path for 2-pin nets (common case): avoid generic sort/unique.
        if (coarse_nodes_buf.size() == 2) {
            NodeID a = coarse_nodes_buf[0];
            NodeID b = coarse_nodes_buf[1];
            if (a == b) {
                continue;
            }
            if (a > b) {
                std::swap(a, b);
            }
            coarse_nodes_buf[0] = a;
            coarse_nodes_buf[1] = b;
        } else {
            // Sort and remove duplicates (faster than unordered_set for small sizes)
            std::sort(coarse_nodes_buf.begin(), coarse_nodes_buf.end());
            auto last = std::unique(coarse_nodes_buf.begin(), coarse_nodes_buf.end());
            coarse_nodes_buf.erase(last, coarse_nodes_buf.end());
        }
        if (coarse_nodes_buf.size() <= 1) continue;

        // Check for parallel net by hash bucket + exact compare.
        const size_t key_hash = vector_hash(coarse_nodes_buf);
        auto bucket_it = hash_to_indices.find(key_hash);
        size_t matched_idx = static_cast<size_t>(-1);
        if (bucket_it != hash_to_indices.end()) {
            const std::vector<size_t>& candidates = bucket_it->second;
            for (size_t idx : candidates) {
                const std::vector<NodeID>& existing_key = net_keys[idx];
                if (existing_key.size() == coarse_nodes_buf.size() &&
                    std::equal(existing_key.begin(), existing_key.end(),
                               coarse_nodes_buf.begin())) {
                    matched_idx = idx;
                    break;
                }
            }
        }

        if (matched_idx != static_cast<size_t>(-1)) {
            // Parallel net found, merge weights.
            net_weights[matched_idx] += fine_hg.getNetWeight(fine_net);
        } else {
            // New net (preserve first-seen order for deterministic net IDs).
            const size_t idx = net_keys.size();
            net_keys.push_back(coarse_nodes_buf);
            net_weights.push_back(fine_hg.getNetWeight(fine_net));
            auto insert_result = hash_to_indices.emplace(key_hash, std::vector<size_t>());
            insert_result.first->second.push_back(idx);
        }
    }
    
    // Add nets to coarse graph
    coarse_hg.reserveNets(net_keys.size());
    for (size_t i = 0; i < net_keys.size(); ++i) {
        EdgeID coarse_net = coarse_hg.addNet(net_weights[i], false);
        for (NodeID coarse_node : net_keys[i]) {
            coarse_hg.addNodeToNet(coarse_net, coarse_node);
        }
    }
}

bool ClusterMatching::adjMatchAreaJudging(const Hypergraph& hg, NodeID node1, NodeID node2, Weight max_weight) const {
    Weight node1_weight = hg.getNodeWeight(node1);
    Weight node2_weight = hg.getNodeWeight(node2);
    if (node1_weight + node2_weight > max_weight) {
        return false;
    }
    return true;
}
} // namespace consmlp
