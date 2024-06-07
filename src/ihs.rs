use crate::graph;
use crate::{sifting,heuristic,quickhs};
use crate::scc::SCC;
use std::collections::{HashSet,HashMap};
use itertools::{ Itertools, iproduct };
use bit_set::BitSet;
use rand::Rng;
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::time::Instant;
use highs::RowProblem;

/// Returns an ordering of the vertices such that the directed edges x -> y violating this ordering (x comes after y) form a *minimum* weight [feedback arc set] of the given graph 
/// (i.e. a topological ordering of the spanning DAG of the FAS), which is usually a strongly connected component of the input.
///
/// This function implements the [implicit hitting set] algorithm: We consider a *dynamic* set system, in which the sets are the cycles
/// of the input graph. Since there may be exponentially many such cycles, we start with a small initial subset of them. We find a
/// hitting set of these cycles and check whether it happens to hit *all* cycles of the input graph. If so, we found an optimal feedback arc set,
/// since we computed a minimum hitting set. If not, we search for cycles in the remaining graph, add them to our set system, and recompute
/// the hitting set.
///
/// # Computation of the Hitting Set
/// The hitting set is computed in each round with `hitting_set`, which currently uses the LP/ILP solver [highs].
///
/// [feedback arc set]: https://en.wikipedia.org/wiki/Feedback_arc_set
/// [implicit hitting set]: https://link.springer.com/chapter/10.1007/978-3-642-13509-5_14
/// [highs]: https://highs.dev/
pub fn solve(scc: &SCC) -> Vec<usize> {  
    let verbose = false; 
    let start = Instant::now();
    let     n = scc.n; // number of nodes  
    let mut m = 0;       // number of arcs, arc ids are 0...m-1
    let mut rng = thread_rng();

    // computing arc ids (i.e., variables) and their weight
    let mut weights = Vec::new();
    let mut idarc = Vec::new(); 
    let mut arcs    = HashMap::<(usize,usize), usize>::new();
    iproduct!(0..n, 0..n)
	.filter(|(u,v)| scc.w[*u][*v] > 0 )
	.for_each(|(u, v)| { weights.push(scc.w[u][v]); arcs.insert((u,v), m); idarc.push((u, v)); m += 1; } );
    let arcid = |u: usize, v: usize| { *arcs.get(&(u,v)).unwrap() };
    
    // map cycles (as set of vertices) to a (sorted) set of arc ids
    let cycle2arcs = |cycle: Vec<usize>| {
	let mut c = cycle.into_iter().circular_tuple_windows::<(usize,usize)>().map(|uv| arcid(uv.0,uv.1)).collect_vec();
	c.sort(); c	
    };

    let mut cycles;
    if scc.n >= 1000 {
        cycles = graph::three_cycles(&scc.g).into_iter().map(cycle2arcs).collect_vec();
    } else {
        cycles = init_cycles_sifting(scc).into_iter().map(cycle2arcs).collect_vec();
    }

    // constrain number of initial cycles added
    let init_cap = 100_000;
    if cycles.len() > init_cap {
        cycles.shuffle(&mut rng);
        cycles.truncate(init_cap);
    }
    if verbose {
        println!("#cols {}", weights.len());
        println!("#rows init {}", cycles.len());
    }

    // repeatedly find a feedback arc set (a hitting set on the cycles) and compute new (still not hit) cycles
    let mut fas;
    let mut sol;
    let mut vals;

    loop {
	// This (outer) loop is the main ILP loop. If it breaks, an optimal fas was found.
	loop {
	    // This (inner) loop finds hitting sets using an LP relaxation.
	    // If it breaks, we have a hitting set for all cycles, but it is not necessarily optimal. // This print statement can be helpful for debugging.
            // println!("#rows {}", cycles.len());
            loop {
		// Quick heuristic loop.
		// This loop computes a hitting set using a simple max-degree heuristic.
		// If it breaks, we have a hitting set for all cycles, but it is not necessarily optimal.
                if verbose {
                    println!("#rows heur {} elapsed {}", cycles.len(), start.elapsed().as_secs());
                }
		let fas = quickhs::hitting_set(&cycles, &weights);
                let oldsz = cycles.len();
                new_cycles(&scc.w, &scc.g, &fas, arcid).into_iter()
                    .map(cycle2arcs)
                    .for_each(|c| cycles.push(c) );
                if cycles.len() == oldsz {
                    break;
                }
	    }	  
	    // This (inner) loop finds hitting sets using an LP relaxation.
	    // If it breaks, we have a hitting set for all cycles, but it is not necessarily optimal.
            if verbose {
                println!("#rows lp {} elapsed {}", cycles.len(), start.elapsed().as_secs());
            }
	    (fas, sol) = hitting_set(&cycles, &weights, true);
	    vals = sol.columns().to_vec(); 
            
            if verbose {
                let mut lb = 0.0;
                let mut ub = 0.0;
                for i in 0..vals.len() {
                    lb += vals[i] * weights[i] as f64;
                    if fas.contains(i) {
                        ub += weights[i] as f64;
                    }
                }
                println!("lb {} ub {} elapsed {}", lb, ub, start.elapsed().as_secs());
            }

            let oldsz = cycles.len();
            new_cycles(&scc.w, &scc.g, &fas, arcid).into_iter()
                .map(cycle2arcs)
                .for_each(|c| cycles.push(c) );
	    if cycles.len() == oldsz {
                break;
            }
	}
        // compute bounds based on last LP solution and rounding
        let mut lb = 0.0;
        let mut ub = 0.0;
        for i in 0..vals.len() {
            lb += vals[i] * weights[i] as f64;
            if fas.contains(i) {
                ub += weights[i] as f64;
            }
        }
        // if bounds match by less than full integer, the rounded LP relaxation is optimal
        if ub - lb < 0.99 {
            let resfas = (0..n).flat_map( |v| std::iter::once(v).cartesian_product(scc.g[v].iter()) )
                .filter    ( |(u,v)| fas.contains(arcid(*u,**v))                                    )
                .map(|(u,v)| (u,*v)).collect_vec();
            return scc.fas_to_ordering(&resfas);
        }

        // try to add some extra cycles if ILP would be started
        let oldsz = cycles.len();
        let mut extra_cycles = HashSet::new();
        for idx in 0..m {
            let (j, i) = idarc[idx];
            if vals[idx] > 0.01 && vals[idx] < 0.99 {
                graph::find_more_paths_relaxed(&scc.g, i, j, 10, &vals, arcid, cycle2arcs).into_iter()
                .map(cycle2arcs)
                .for_each(|c| {extra_cycles.insert(c);} );
            }
        }
        for c in extra_cycles.into_iter() {
            cycles.push(c);
        }
        if cycles.len() > oldsz {
            continue;
        }
        
        if verbose { 
            println!("start ilp solve {}", start.elapsed().as_secs());
        }
        // else, solve the ILP
	(fas, _) = hitting_set(&cycles, &weights, false);
        let oldsz = cycles.len();
        new_cycles(&scc.w, &scc.g, &fas, arcid).into_iter()
            .map(cycle2arcs)
            .for_each(|c| cycles.push(c));
        if cycles.len() == oldsz {
            break;
        }
    }
    
    let resfas = (0..n).flat_map( |v| std::iter::once(v).cartesian_product(scc.g[v].iter()) )
	.filter    ( |(u,v)| fas.contains(arcid(*u,**v))                                    )
	.map(|(u,v)| (u,*v)).collect_vec();
    scc.fas_to_ordering(&resfas)
}

pub fn generate_new_cycles(fas: &Vec<(usize, usize)>, h: &Vec<Vec<usize>>) -> Vec<Vec<usize>> { 
    let sccs = graph::compute_sccs(h);
    let mut to_scc = vec![0; h.len()];
    for i in 0..sccs.len() {
        for v in sccs[i].iter().cloned() {
            to_scc[v] = i;
        }
    }
    let mut new_cycles: HashSet<Vec<usize>> = HashSet::new();
    for (u,v) in fas.iter().cloned() {
        if to_scc[u] != to_scc[v] {
            continue;
        }
        let k = 10;
        graph::find_more_cycles_bfs(h, &to_scc, v, u, k).into_iter().for_each(|c| {new_cycles.insert(c);}); 
    }
    new_cycles.into_iter().collect()
}

pub fn find_fas(g: &Vec<Vec<usize>>, w: &Vec<Vec<u64>>) -> Vec<(usize, usize)> {
    let mut fas: Vec<(usize, usize)> = Vec::new();
    let sccs = graph::compute_sccs(g);
    for scc in sccs.iter().cloned() {
        if scc.len() == 1 {
            continue;
        }
        // construct SCC 
        let mut scc_w = vec![vec![0; scc.len()]; scc.len()];
        for i in 0..scc.len() {
            for j in 0..scc.len() {
                scc_w[i][j] = w[scc[i]][scc[j]];
            }   
        }
        let scc_g = graph::get_subgraph(g, &scc);
        let scc_obj = SCC::new(scc, scc_w, scc_g);
        let scc_ordering = heuristic::kahns_heuristic(&scc_obj);
        let scc_fas = scc_obj.ordering_to_fas(&scc_ordering);
        for (u, v) in scc_fas {
            fas.push((scc_obj.labels[u], scc_obj.labels[v]));
        }
    }
    fas
}

fn init_cycles_sifting(scc: &SCC) -> Vec<Vec<usize>> {
    let mut new_cycles: HashSet<Vec<usize>> = HashSet::new();
    let fas = find_fas_initial(scc); 
    let to_scc = vec![0; scc.g.len()];
    for (u,v) in fas.iter().cloned() {
        if to_scc[u] != to_scc[v] {
            continue;
        }
        let k = 20; 
        graph::find_k_cycles_bfs(&scc.g, &to_scc, v, u, k).into_iter().for_each(|c| {new_cycles.insert(c);}); 
    }
    new_cycles.into_iter().collect()
}

// only for initial cycle generation
pub fn find_fas_initial(scc: &SCC) -> Vec<(usize, usize)> {
    let sccs = vec![SCC::new(scc.labels.clone(), scc.w.clone(), scc.g.clone())]; 
    let ordering = sifting::hillclimber_sifting(&sccs, sifting::insertionplus_sifting(&sccs));
    sccs[0].ordering_to_fas(&ordering[0])
}

/// Computes a set of cycles that was not hit by the given feedback arc set.
fn new_cycles(w: &Vec<Vec<u64>>, g: &Vec<Vec<usize>>, fas: &BitSet, arcid: impl Fn(usize,usize)->usize) -> Vec<Vec<usize>> {
    let mut residual: Vec<Vec<usize>> = vec![Vec::new(); g.len()];
    (0..w.len()).flat_map( |v| std::iter::once(v).cartesian_product(g[v].iter()) )
	.filter    ( |(u,v)| !fas.contains(arcid(*u,**v))                        )
	.for_each  ( |(u,v)|  residual[u].push(*v)                               );

    let mut cycles = generate_new_cycles( 
	&find_fas(&residual, w), 
	&residual
    );
  
    let update_cap = 20_000;
    if cycles.len() > update_cap {
        cycles.shuffle(&mut thread_rng());
        cycles.truncate(update_cap);
    }
    cycles
}

/// Computes an minimum weight hitting set of the given set system.
///
/// The sets should be ordered and contain elements from `0...weights.len()-1`.
/// Correspondingly, the weight of element `i` is `weights[i]`.
///
/// This function finds the hitting set with the trivial ILP formulation using the highs lp solver.
fn hitting_set(sets: &Vec<Vec<usize>>, weights: &Vec<u64>, relaxation: bool) -> (BitSet, highs::Solution) {
    let mut problem = RowProblem::new();
    let mut cols = Vec::new();
    for weight in weights {
	match relaxation {
	    false => cols.push( problem.add_integer_column(*weight as f64, 0..=1) ),
	    true  => cols.push( problem.add_column(        *weight as f64, 0..=1) ),

	};
    }

    for set in sets.iter() {
        problem.add_row(1..,  &set.iter().map(|i| (cols[*i], 1.0)).collect_vec())
    }

    let mut model = problem.optimise(highs::Sense::Minimise);

    match relaxation {
        false => model.set_option("solver", "choose"),
        true =>  model.set_option("solver", "ipm") 
    };
    
    model.set_option("parallel", "off");
    model.set_option("threads", 1);

    let solved_model = model.solve();
    let sol = solved_model.get_solution();
  
    let hs;
    let cols = sol.columns().to_vec();
    if relaxation { 
        let mut rng = thread_rng();
        // first add elements to hs which have value 1.0
        let mut safe_hs = BitSet::new();
        let mut safe_hit = vec![false; sets.len()];
        for i in 0..sets.len() {
            for j in sets[i].iter().cloned() {
                if cols[j] > 0.99 {
                    safe_hs.insert(j); // should not insert duplicates
                    safe_hit[i] = true;
                    break;
                }
            }
        }
        let mut new_sets = Vec::new();
        for i in 0..sets.len() {
            if !safe_hit[i] {
                new_sets.push(sets[i].clone());
            }
        }
        if new_sets.is_empty() {
            hs = safe_hs;
        } else {
            // construct remaining instance
            let mut memb_graph: Vec<Vec<usize>> = vec![Vec::new(); weights.len()];
            for i in 0..new_sets.len() {
                for j in new_sets[i].iter().cloned() {
                    memb_graph[j].push(i);
                }
            }
            let mut rest_hs: Vec<usize> = Vec::new(); 
            for set in new_sets.iter() {
                if set.iter().any(|e| rest_hs.contains(e)) { continue; }
                if let Some((_,e)) = set.iter().map(|e| (cols[*e] + rng.gen::<f64>() / 1000.0,e)).max_by(|a, b| a.0.total_cmp(&b.0)) {
                    rest_hs.push(*e);
                }
            }
            let mut hit_by: Vec<Vec<usize>> = vec![Vec::new(); new_sets.len()];
            for i in rest_hs.iter().cloned() {
                for j in memb_graph[i].iter().cloned() {
                    hit_by[j].push(i);
                }
            } 
            let mut curval: i64 = 0;
            let mut curhs: Vec<usize> = rest_hs.clone();
            let mut optval: i64 = curval;
            let mut opths: Vec<usize> = rest_hs.clone();
            
            let temperature: f64 = 1.0;
            let iterations = std::cmp::min(new_sets.len() * 1_000, 1_000_000);
            for _ in 0..iterations {
                let mut kick_pos;
                let mut kick_out;
                loop {
                    kick_pos = rng.gen_range(0..curhs.len());
                    kick_out = curhs[kick_pos];
                    if rng.gen::<f64>() < 1.0 - cols[kick_out] {
                        break;
                    }
                }
                let mut delta = -(weights[kick_out] as i64);
                curhs.swap_remove(kick_pos);
                let mut put_in: Vec<usize> = Vec::new();
                for i in memb_graph[kick_out].iter().cloned() {
                    hit_by[i].retain(|x| *x != kick_out); 
                    if hit_by[i].len() == 0 {
                        let mut mx = 0.0;
                        let mut mxi = 0;
                        for j in new_sets[i].iter().cloned() {
                            if j == kick_out {
                                continue;
                            } 
                            if mx < cols[j] + rng.gen::<f64>() / 1000.0 {
                                mx = cols[j];
                                mxi = j;
                            }
                        }
                        delta += weights[mxi] as i64;
                        curhs.push(mxi);
                        put_in.push(mxi);
                        for k in memb_graph[mxi].iter().cloned() {
                            hit_by[k].push(mxi);
                        }
                    }
                }
                if (-delta as f64 / temperature).exp() >= rng.gen() {
                    curval += delta;
                    if curval < optval {
                        opths = curhs.clone();
                        optval = curval;
                    }
                } else {
                    curhs.retain(|x| !put_in.contains(x));
                    for i in put_in.iter().cloned() {
                        for k in memb_graph[i].iter().cloned() {
                            hit_by[k].retain(|x| *x != i);
                        }
                    }
                    curhs.push(kick_out);
                    for k in memb_graph[kick_out].iter().cloned() {
                        hit_by[k].push(kick_out);
                    }
                }
            }
            for i in opths.iter().cloned() {
                safe_hs.insert(i);
            }
            hs = safe_hs;
        }
    } else {
        let mut round_hs = BitSet::new();
	(0..weights.len()).filter(|i| cols[*i].round() as i64 == 1)
	    .for_each(|i| { round_hs.insert(i); } );
        hs = round_hs;
    }
    (hs, sol)   
}

