use crate::bipartite_graph::BipartiteGraph;
use crate::ihs;

/// Entry point for solvers based on the reduction to weighted directed feedback arc set (WDFAS).
/// Performs the reduction and calls the solver for each SCC of the obtained instance. 
/// Aftewards the results of the individual calls are glued together for the final result.
pub fn wdfas(g: &BipartiteGraph) -> Vec<usize> {
    let sccs = g.reduce();
    let mut ordering = Vec::new();
    for u in g.isolated.iter().cloned() {
        ordering.push(u);
    }
    for scc in sccs.iter() {
        if scc.n == 1 {
            for twin in g.ids[scc.labels[0]].iter().cloned() {
                ordering.push(twin);
            }
        } else {
            let scc_ordering = ihs::solve(scc);
            for v in scc_ordering.iter().cloned() {
                for twin in g.ids[scc.labels[v]].iter().cloned() {
                    ordering.push(twin);
                }
            }
        }
    }
    for el in &mut ordering {
        *el += g.n0 + 1;
    }
    ordering
}
