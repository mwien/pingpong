use crate::graph;

pub struct SCC {
    pub n: usize,
    pub labels: Vec<usize>,
    pub w: Vec<Vec<u64>>,
    pub g: Vec<Vec<usize>> 
}

impl SCC {
    pub fn new(labels: Vec<usize>, w: Vec<Vec<u64>>, g: Vec<Vec<usize>>) -> SCC {
        SCC {n: labels.len(), labels, w, g }
    }
    pub fn fas_to_ordering(&self, fas: &Vec<(usize, usize)>) -> Vec<usize> {
        let mut fas_lookup: Vec<Vec<bool>> = vec![vec![false; self.n]; self.n]; 
        for (u,v) in fas.iter() {
            fas_lookup[*u][*v] = true;
        }
        // construct rest graph with edges not in feedback arc set 
        // start with subgraph of h induced by scc
        let mut rest_graph: Vec<Vec<usize>> = vec![Vec::new(); self.n];
        for u in 0..self.n { 
            for v in self.g[u].iter().cloned() {
                if !fas_lookup[u][v] { 
                    rest_graph[u].push(v);
                }
            }
        }
        graph::top_ordering(&rest_graph)
    }
    
    pub fn ordering_to_fas(&self, ordering: &Vec<usize>) -> Vec<(usize, usize)> { 
        let mut fas = Vec::new();
        let mut invorder = vec![0; ordering.len()];
        for i in 0..ordering.len() {
            invorder[ordering[i]] = i;
        }
        for u in 0..self.n {
            for v in self.g[u].iter().cloned() {
                if invorder[v] < invorder[u] {
                    fas.push((u, v));
                }
            }
        }
        fas 
    }
}
