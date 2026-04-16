#ifndef ZXPART_MULTILEVEL_COARSENER_H
#define ZXPART_MULTILEVEL_COARSENER_H

#include "Coarsener.h"
#include "datastructures/HypergraphHierarchy.h"
#include <memory>
#include <vector>

namespace consmlp {

/**
 * @brief Manages multilevel coarsening process
 * 
 * Features:
 * - Iterative coarsening until threshold
 * - Contraction ratio control (1.5-1.7)
 * - Statistics collection
 */
class MultilevelCoarsener {
public:
    /**
     * @brief Constructor
     * @param coarsener Coarsening algorithm to use
     */
    explicit MultilevelCoarsener(std::unique_ptr<Coarsener> coarsener);
    
    /**
     * @brief Perform multilevel coarsening
     * @param finest_hg Finest hypergraph
     * @param silent If true, suppress progress output
     * @return Hypergraph hierarchy
     */
    HypergraphHierarchy coarsen(Hypergraph&& finest_hg, bool silent = false);
    
    /**
     * @brief Get coarsening statistics for all levels
     * @return Vector of statistics
     */
    const std::vector<CoarseningStats>& getStatistics() const {
        return stats_;
    }
    
private:
    std::unique_ptr<Coarsener> coarsener_;
    std::vector<CoarseningStats> stats_;
};

} // namespace consmlp

#endif // ZXPART_MULTILEVEL_COARSENER_H
