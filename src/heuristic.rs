use crate::scc::SCC;

pub fn kahns_heuristic(scc: &SCC) -> Vec<usize> {
    let mut inweightsum: Vec<u64> = vec![0; scc.labels.len()];
    for u in 0..scc.labels.len() {
        for v in scc.g[u].iter().cloned() {
            inweightsum[v] += scc.w[u][v];
        }
    } 
    let mut order: Vec<usize> = vec![0; scc.labels.len()];
    for i in 0..scc.labels.len() {
        let mut mn: f64 = f64::MAX;
        let mut mnj: usize = 0;
        for j in 0..scc.labels.len() {
            let val = inweightsum[j] as f64;
            if val < mn {
                mn = val;
                mnj = j;
            }
        }
        order[i] = mnj;
        inweightsum[mnj] = u64::MAX; // enforce mnj is not chosen again
        for j in scc.g[mnj].iter().cloned() {
            inweightsum[j] -= scc.w[mnj][j];
        }
    }
    order 
}
