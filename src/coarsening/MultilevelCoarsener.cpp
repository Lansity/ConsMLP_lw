#include "coarsening/MultilevelCoarsener.h"
#include <iostream>

namespace consmlp {

MultilevelCoarsener::MultilevelCoarsener(std::unique_ptr<Coarsener> coarsener)
    : coarsener_(std::move(coarsener))
{
}

HypergraphHierarchy MultilevelCoarsener::coarsen(Hypergraph&& finest_hg, bool silent) {
    stats_.clear();
    
    // Create hierarchy with finest level
    HypergraphHierarchy hierarchy(std::move(finest_hg));
    
    uint32_t level = 0;
    
    if (!silent) {
        std::cout << "Starting multilevel coarsening..." << std::endl;
    }
    
    // Iteratively coarsen
    while (true) {
        const Hypergraph& current_hg = hierarchy.getLevel(level).getHypergraph();
        NodeID num_nodes = current_hg.getNumNodes();
        
        // Check stopping criterion
        if (coarsener_->shouldStopCoarsening(num_nodes)) {
            break;
        }
        
        // Perform coarsening
        CoarseningStats stats = coarsener_->coarsen(hierarchy, level);
        stats_.push_back(stats);
        
        if (stats.coarse_nodes < 1000) {
            if (stats.coarse_nodes < 100 || stats.coarse_nets < 200) {
                break;
            }
            int valued_nodes = 0;
            auto coarse_hg = hierarchy.getLevel(level + 1).getHypergraph();
            for (NodeID node = 0; node < stats.coarse_nodes; ++node) {
                if (coarse_hg.getNodeWeight(node) > 0 && coarse_hg.getNodeType(node) != NodeType::IO && !coarse_hg.isNodeIgnored(node) && !coarse_hg.isNodeFixed(node)) {
                    valued_nodes++;
                }
            }
            if (valued_nodes < 80) {
                break;
            }
        }
        if (stats.contraction_ratio < 1.1) {
            break;
        }
        
        level++;
        
        // Safety limit on number of levels
        if (level >= 40) {
            std::cout << "Stopping: Maximum levels reached" << std::endl;
            break;
        }
    }
    
    return hierarchy;
}

} // namespace consmlp
