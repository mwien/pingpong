use std::io::BufRead;
use std::collections::HashMap;
use std::error::Error;
use crate::graph;
use crate::scc::SCC;
use std::cmp;

pub struct BipartiteGraph {
    pub n0: usize,
    pub n1: usize, 
    pub adjs: Vec<Vec<usize>>, // neighbors of vertices 1, ..., n1 in second partition
    pub ids: Vec<Vec<usize>>, // original ids of vertices 1, ..., n1
    pub isolated: Vec<usize>, // ids of additional isolated vertices
}

impl BipartiteGraph {
    pub fn new(n0: usize, n1: usize, adjs: Vec<Vec<usize>>, ids: Vec<Vec<usize>>, isolated: Vec<usize>) -> BipartiteGraph {
        BipartiteGraph { n0, n1, adjs, ids, isolated }
    }

    pub fn new_from_stdin() -> Result<BipartiteGraph, Box<dyn Error>> {
        let mut ingraph: Option<Vec<Vec<usize>>> = None;
        let mut n0: usize = 0;
        for line in std::io::stdin().lock().lines() {
            let line = line?;
            let ll: Vec<&str> = line.split(' ').collect();
            match ll[0] {
                "c" => {} // skip comments 
                "p" => { // parse header 
                    let a = ll[2].parse::<usize>()?;
                    let b = ll[3].parse::<usize>()?;
                    n0 = a;
                    ingraph = Some(vec![Vec::new(); b]);
                },
                _ => { // parse                    
                    match ingraph {
                        None => return Err(From::from("c Found edge before p-line. Abort!")),
                        Some(ref mut ingraph) => {
                            if ll.len() < 2 {
                                continue;
                            }
                            let a = ll[0].parse::<usize>()?;
                            let b = ll[1].parse::<usize>()?;
                            ingraph[b-n0-1].push(a-1);
                        }
                    }
                }
            }
        }
        match ingraph {
            Some(mut ingraph) => {
                let mut isolated: Vec<usize> = Vec::new();
                let mut adjtotwins: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
                for i in 0..ingraph.len() {
                    if ingraph[i].is_empty() {
                        isolated.push(i);
                    } else {
                        ingraph[i].sort();
                        adjtotwins.entry(ingraph[i].clone()).and_modify(|twins| twins.push(i)).or_insert(vec![i]);
                    }
                }
                let (adjs, ids): (Vec<Vec<usize>>, Vec<Vec<usize>>) = adjtotwins.into_iter().unzip();
                Ok(BipartiteGraph::new(n0, ids.len(), adjs, ids, isolated))
            },
            None => Err(From::from("c Failed to parse a graph! Maybe the input was empty?"))
        }
    }

    #[inline(always)]
    pub fn pair_crossing_number(&self, u: usize, v: usize) -> u64 {
        let nu = &self.adjs[u];
        let nv = &self.adjs[v];
        let mut cn: u64 = 0;
        if nu.len() <= nv.len() {
            for x in nu.iter() {
                cn += match nv.binary_search(x) {
                    Ok(i) =>  i as u64,
                    Err(i) => i as u64
                };
            }
        } else {
            for x in nv.iter() {
                cn += nu.len() as u64 - match nu.binary_search(x) {
                    Ok(i) =>  i as u64 + 1, 
                    Err(i) => i as u64
                }; 
            }
        }
        cn
    }

    pub fn crossing_matrix(&self) -> Vec<Vec<u64>> {
        let mut cm = vec![vec![0; self.n1]; self.n1];
        for u in 0..self.n1 {
            for v in 0..self.n1 {
                if u == v {
                    continue;
                }
                cm[u][v] = self.pair_crossing_number(u, v) * self.ids[u].len() as u64 * self.ids[v].len() as u64;
            }
        }
        cm
    }

    pub fn reduce(&self) -> Vec<SCC> {
        let mut cm = self.crossing_matrix();
        for u in 0..cm.len() {
            for v in u+1..cm.len() {
                let mn: u64 = cmp::min(cm[u][v], cm[v][u]);
                cm[u][v] -= mn;
                cm[v][u] -= mn;
            }
        } 
        let mut h: Vec<Vec<usize>> = vec![Vec::new(); cm.len()];
        for u in 0..cm.len() {
            for v in 0..cm.len() {
                if cm[v][u] != 0 {
                    h[u].push(v); 
                }
            }
        } 
        let sccs = graph::compute_sccs(&h);
        let mut result: Vec<SCC> = Vec::new(); 
        for scc in sccs.iter().cloned() {
            let mut w: Vec<Vec<u64>> = vec![vec![0; scc.len()]; scc.len()];
            let mut g: Vec<Vec<usize>> = vec![Vec::new(); scc.len()];
            for i in 0..scc.len() {
                for j in 0..scc.len() {
                    w[i][j] = cm[scc[j]][scc[i]];
                    if w[i][j] != 0 {
                        g[i].push(j);
                    }
                }
            }
            result.push(SCC::new(scc, w, g));
        }
        result
    }
}
