use bit_set::BitSet;
use std::collections::HashSet;

pub fn hitting_set(sets: &Vec<Vec<usize>>, weights: &Vec<u64>) -> BitSet {
    // The degree of an element is just the number of sets it is in.
    let mut watch = vec![HashSet::<usize>::default(); weights.len()];
    for (i,set) in sets.iter().enumerate() {
	set.iter().for_each(|e| { watch[*e].insert(i); } );
    }

    // A bucket priority queue.
    // For every possible degree, we have a bucket of all elements of that degree.
    // Degrees can only be reduced.
    let mut queue_head   = watch.iter().map(|w| w.len()).max().unwrap();
    let mut bucket_queue = vec![Vec::new(); queue_head+1];
    for e in 0..weights.len() {
	bucket_queue[ watch[e].len() ].push(e);
    }
    
    // As long as there are elements in sets that are not hit, take one of them.
    let mut hs = BitSet::new();
    while queue_head > 0 {
	if let Some(e) = bucket_queue[ queue_head ].pop() {	
	    if watch[e].len() != queue_head { continue; } // degree of e has changed
	    hs.insert(e);
	    let ids = watch[e].iter().cloned().collect::<Vec<usize>>();
	    for set_id in ids {
		for v in sets[set_id].iter() {
		    watch[*v].remove(&set_id);
		    bucket_queue[ watch[*v].len() ].push(*v);
		}
	    }
	} else { queue_head -= 1; }
    }
    hs
}
