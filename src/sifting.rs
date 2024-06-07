use crate::scc::SCC;
use rand::thread_rng;
use rand::seq::SliceRandom;
use rand::Rng;

fn get_inv_w(sccs: &Vec<SCC>) -> Vec<Vec<Vec<u64>>> {
    let mut inv_w: Vec<Vec<Vec<u64>>> = Vec::new();
    for scc in sccs {
        let mut scc_inv_w = vec![vec![0; scc.n]; scc.n];
        for j in 0..scc.n {
            for k in 0..scc.n {
                scc_inv_w[j][k] = scc.w[k][j];
            }
        }
        inv_w.push(scc_inv_w);
    }
    inv_w 
}

fn insert_cost_per_pos(scc: &SCC, diffs: &Vec<Vec<u64>>, perm: &Vec<usize>, v: usize) -> Vec<u64> {
    let n = perm.len();
    let mut pre: Vec<u64> = vec![0; n+1];
    let mut suf: Vec<u64> = vec![0; n+1];
    for i in 0..n {
        pre[i+1] = pre[i] + scc.w[v][perm[i]];
    }
    let mut perm_rev = perm.clone();
    perm_rev.reverse();
    for i in 0..n {
        suf[i+1] = suf[i] + diffs[v][perm_rev[i]];
    }
    suf.reverse();
    pre.iter().zip(&suf).map(|(p, s)| p + s).collect()
}

fn get_min(val: &Vec<u64>) -> (u64, Vec<usize>) {
    let minval = *val.iter().min().unwrap();
    let minima = val.iter()
        .enumerate()
        .filter(|(_, &x)| x == minval)
        .map(|(idx, _)| idx)
        .collect();
    (minval, minima)
}

pub fn insertionplus_sifting(sccs: &Vec<SCC>) -> Vec<Vec<usize>> {
    let inv_w = get_inv_w(sccs); // could precompute this
    let mut ordering: Vec<Vec<usize>> = Vec::new();
    let mut rng = thread_rng();
    for i in 0..sccs.len() {
        let scc = &sccs[i];
        let mut scc_ordering: Vec<usize> = Vec::new();
        let mut vertices: Vec<usize> = (0..scc.n).collect();
        vertices.shuffle(&mut rng);
        for v in vertices.iter().cloned() {
            let cost = insert_cost_per_pos(scc, &inv_w[i], &mut scc_ordering, v);
            let (_, minima) = get_min(&cost);
            scc_ordering.insert(*minima.choose(&mut rng).unwrap(), v);
            if scc_ordering.len() % 50 == 0 {
                // put into function
                let mut iter = 0;
                let mut last_improvement = 0;
                while iter - last_improvement < 2*scc_ordering.len() { 
                    let vpos = rng.gen_range(0..scc_ordering.len());
                    let v = scc_ordering[vpos];
                    scc_ordering.remove(vpos); 
                    let cost = insert_cost_per_pos(scc, &inv_w[i], &scc_ordering, v);
                    let (mincost, minima) = get_min(&cost);
                    let previous_cost = cost[vpos]; 
                    let delta = previous_cost - mincost;
                    if delta > 0 {
                        last_improvement = iter;
                    }
                    let mut inspos;
                    loop {
                        inspos = *minima.choose(&mut rng).unwrap();
                        if inspos !=  vpos || minima.len() == 1 {
                            break;
                        }
                    }
                    scc_ordering.insert(inspos, v);
                    iter += 1;
                }
                
            }
        }
        ordering.push(scc_ordering);
    }
    ordering
    
}

pub fn hillclimber_sifting(sccs: &Vec<SCC>, initial_ordering: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let inv_w = get_inv_w(sccs); // could precompute this
    let mut ordering = initial_ordering;
    let mut rng = thread_rng();
    let mut iter = 0;
    let mut last_improvement = 0;
    while iter - last_improvement < 4 { 
        for i in 0..sccs.len() {
            let scc = &sccs[i];
            if scc.n == 1 { continue; }
            let scc_ordering = &mut ordering[i];            
            // do shuffles or just take random elements?
            let mut vertices: Vec<usize> = (0..scc.n).collect();
            vertices.shuffle(&mut rng);
            for v in vertices.iter().cloned() {
                // for now always remove element first, later optimize this
                let vpos = scc_ordering.iter().position(|&x| x == v).unwrap();
                scc_ordering.remove(vpos); 
                let cost = insert_cost_per_pos(scc, &inv_w[i], scc_ordering, v);
                let (mincost, minima) = get_min(&cost);
                let previous_cost = cost[vpos]; 
                let delta = previous_cost - mincost;
                if delta > 0 {
                    last_improvement = iter;
                }
                let mut inspos;
                loop {
                    inspos = *minima.choose(&mut rng).unwrap();
                    if inspos !=  vpos || minima.len() == 1 {
                        break;
                    }
                }
                scc_ordering.insert(inspos, v);
            }
        }
        iter += 1;
    }
    ordering
}
