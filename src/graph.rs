use std::collections::{VecDeque, HashSet};
use rand::seq::SliceRandom;
use rand::thread_rng;

// This crate contains pure graph functionality.
// Graphs are represented as Vec<Vec<usize>>. 

/// Constructs subgraph of g induced by subset. In the resulting graph vertex i corresponds to
/// vertex subset[i] in the original graph.
pub fn get_subgraph(g: &Vec<Vec<usize>>, subset: &Vec<usize>) -> Vec<Vec<usize>> {
    let mut imp: Vec<i32> = vec![-1; g.len()]; 
    for i in 0..subset.len() {
        imp[subset[i]] = i as i32;
    }
    let mut h: Vec<Vec<usize>> = vec![Vec::new(); subset.len()];
    for i in 0..subset.len() {
        for v in g[subset[i]].iter().cloned() {
            let newv = imp[v];
            if newv != -1 {
                h[i].push(newv as usize);
            }
        }
    }
    h
}

fn top_ordering_dfs(g: &Vec<Vec<usize>>, vis: &mut Vec<bool>, ord: &mut Vec<usize>, u: usize) {
    if vis[u] {
        return;
    } 
    vis[u] = true;
    for v in g[u].iter().cloned() {
        top_ordering_dfs(g, vis, ord, v);
    }
    ord.push(u);
}

/// Returns topological ordering of directed acyclic graph g. 
pub fn top_ordering(g: &Vec<Vec<usize>>) -> Vec<usize> {
    let mut vis = vec![false; g.len()];
    let mut ord: Vec<usize> = Vec::new();
    for u in 0..g.len() {
        if !vis[u] {
            top_ordering_dfs(g, &mut vis, &mut ord, u);
        }
    }
    ord.reverse();
    ord
}

fn propagate_label(g: &Vec<Vec<usize>>, vertex_labels: &mut Vec<usize>, u: usize, label: usize) {
    if vertex_labels[u] != 0 {
        return;
    }
    vertex_labels[u] = label;
    for v in &g[u] {
        propagate_label(g, vertex_labels, *v, label);
    }
}

/// Returns list of strongly connected components. 
pub fn compute_sccs(h: &Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    // compute top ordering
    let ord = top_ordering(h);
    // get reverse graph
    let mut hrev: Vec<Vec<usize>> = vec![Vec::new(); h.len()];
    for u in 0..h.len() {
        for v in &h[u] {
            hrev[*v].push(u);
        }
    }
    let mut label: usize = 1;
    let mut vertex_labels = vec![0; h.len()];
    for u in ord.into_iter() {
        if vertex_labels[u] == 0 {
            propagate_label(&hrev, &mut vertex_labels, u, label);
            label += 1;
        }
    }
    let mut scc: Vec<Vec<usize>> = vec![Vec::new(); label-1];
    for u in 0..h.len() {
        scc[vertex_labels[u]-1].push(u);    
    }
    scc
}

/// Compute all 3-cycles of g.
pub fn three_cycles(g: &Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let n = g.len();
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    let mut am: Vec<Vec<bool>> = vec![vec![false; n]; n];
    for i in 0..n {
        for j in g[i].iter().cloned() {
            am[i][j] = true;
        }
    }
    for i in 0..n {
        for j in g[i].iter().cloned() {
            if i > j {
                continue;
            }
            for k in g[j].iter().cloned() {
                if i > k {
                    continue;
                }
                if am[k][i] {
                    cycles.push([i, j, k].to_vec());
                }
            }
        }
    }
    cycles
}


fn construct_cycle(pre: &Vec<usize>, i: usize, j: usize) -> Option<Vec<usize>> {
    let mut idx = j;
    let mut cycle: Vec<usize> = Vec::new();
    cycle.push(j);
    while idx != i {
        idx = pre[idx];
        if idx == j {
            return None;
        }
        cycle.push(idx);
    }
    cycle.reverse();
    let mut mn = cycle[0];
    let mut mni: usize = 0;
    for i in 1..cycle.len() {
        if cycle[i] < mn {
            mn = cycle[i];
            mni = i;
        }
    }
    cycle.rotate_left(mni);
    Some(cycle)
}

pub fn find_cycle_bfs(g: &Vec<Vec<usize>>, to_scc: &Vec<usize>, i: usize, j: usize) -> Option<Vec<usize>> {
    let mut rng = thread_rng();
    let mut vis: Vec<bool> = vec![false; g.len()];
    let mut pre: Vec<usize> = vec![0; g.len()];
    let mut q: VecDeque<usize> = VecDeque::new(); 
    vis[i] = true;
    pre[i] = i;
    q.push_back(i);
    while !q.is_empty() {
        let u: usize = q.pop_front().unwrap();
        let mut copyneigh = g[u].clone();
        copyneigh.shuffle(&mut rng);
        for v in copyneigh.iter().cloned() {
            if !vis[v] && to_scc[v] == to_scc[j] {
                q.push_back(v);
                vis[v] = true;
                pre[v] = u;
            }
            if v == j {
                return construct_cycle(&pre, i, j);
            }
        }
    }
    None 
}

pub fn find_k_cycles_bfs(g: &Vec<Vec<usize>>, to_scc: &Vec<usize>, i: usize, j: usize, k: usize) -> Vec<Vec<usize>> { 
    let mut cycles = HashSet::new();
    for _ in 0..k {
        match find_cycle_bfs(g, to_scc, i, j) {
            None => continue,
            Some(c) => cycles.insert(c)
        };
    }
    cycles.into_iter().collect()
}

pub fn find_more_cycles_bfs(g: &Vec<Vec<usize>>, to_scc: &Vec<usize>, i: usize, j: usize, k: usize) -> Vec<Vec<usize>> {
    let mut rng = thread_rng();
    let mut vis: Vec<bool> = vec![false; g.len()];
    let mut pre: Vec<usize> = vec![0; g.len()];
    let mut q: VecDeque<usize> = VecDeque::new(); 
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    vis[i] = true;
    pre[i] = i;
    q.push_back(i);
    while !q.is_empty() {
        let u: usize = q.pop_front().unwrap();
        let mut copyneigh = g[u].clone();
        copyneigh.shuffle(&mut rng);
        for v in copyneigh.iter().cloned() {
            if !vis[v] && to_scc[v] == to_scc[j] {
                q.push_back(v);
                vis[v] = true;
                pre[v] = u;
            }
            if v == j {
                pre[v] = u;
                match construct_cycle(&pre, i, j) {
                    None => (),
                    Some(c) => {
                        cycles.push(c.clone());
                        if cycles.len() > k || c.len() > 4 {
                            return cycles;
                        }
                    }
                };
            }
        }
    }
    cycles
}

pub fn find_more_paths_relaxed(g: &Vec<Vec<usize>>, i: usize, j: usize, k: usize, vals: &Vec<f64>, arcid: impl Fn(usize,usize)->usize, cycle2arcs: impl Fn(Vec<usize>) -> Vec<usize>) -> Vec<Vec<usize>> {
    let mut vis: Vec<bool> = vec![false; g.len()];
    let mut pre: Vec<usize> = vec![0; g.len()];
    let mut q: VecDeque<usize> = VecDeque::new(); 
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    vis[i] = true;
    pre[i] = i;
    q.push_back(i);
    while !q.is_empty() {
        let u: usize = q.pop_front().unwrap();
        // sort by weight 
        let mut indices: Vec<usize> = (0..g[u].len()).collect();
        indices.sort_by(|a, b| vals[arcid(u, g[u][*a])].partial_cmp(&vals[arcid(u, g[u][*b])]).unwrap());
        for idx in indices.iter().cloned() {
            let v = g[u][idx];
            if !vis[v] {
                q.push_back(v);
                vis[v] = true;
                pre[v] = u;
            }
            if v == j {
                pre[v] = u;
                match construct_cycle(&pre, i, j) {
                    None => (),
                    Some(c) => {
                        let mut sm = 0.0;
                        let toarcs = cycle2arcs(c.clone());
                        for ac in toarcs.iter().cloned() {
                            sm += vals[ac];
                        }
                        if sm < 0.99 {
                            cycles.push(c);
                            if cycles.len() > k {
                                return cycles;
                            }
                        }
                    }
                };
            }
        }
    } 
    cycles
}
