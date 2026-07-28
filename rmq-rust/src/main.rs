use std::{
    io::Read,
    path::{Path, PathBuf},
};

fn gauss_summation_interval(l: usize, r: usize) -> usize {
     // This interval is given by (l; r]
     // so the start index is ommited, while the end index
     // is included
    let mut result;
    
    let n = r - l;
    result = n / 2 * (n + 1) + if n & 0b1 == 1 { (n + 1) >> 1 } else { 0 };
    result += l * n;

    result
}

trait Rmq<'a> {
    fn name() -> String;
    /// To save time, only run benchmarks up to this n.
    fn max_n() -> usize {
        usize::MAX
    }
    fn build(data: &'a [u64]) -> Self;
    /// Space usage in bytes.
    fn space(&self) -> usize;
    fn query(&self, l: usize, r: usize) -> u64;
}

/// Trivial implementation that computes each query on the fly.
struct Naive<'a> {
    data: &'a [u64],
}
impl<'a> Rmq<'a> for Naive<'a> {
    fn name() -> String {
        "QuadraticQuery".to_string()
    }
    fn max_n() -> usize {
        // NOTE: Do not use this for the improved implementations!
        10_000
    }
    fn build(data: &'a [u64]) -> Self {
        Self { data }
    }
    fn space(&self) -> usize {
        std::mem::size_of_val(self.data)
    }
    fn query(&self, l: usize, r: usize) -> u64 {
        self.data[l..=r].iter().copied().min().unwrap()
    }
}

// -------------------------------------------------------------
// TODO: Implement the Rmq trait for additional data structures.
// -------------------------------------------------------------

struct LookupTable {
    table: Box<[u64]>,
    array_size: usize,
}
impl<'a> Rmq<'a> for LookupTable {
    fn name() -> String {
        "Lookup table implementation".to_string()
    }

    fn max_n() -> usize {
        // NOTE: Do not use this for the improved implementations!
        10_000
    }
    
    fn build(data: &'a [u64]) -> Self {
        // create array with l * r entries
        // n values for l and n - l for r
        // n / 2 * (n + 1)
        let number_of_entries = gauss_summation_interval(0, data.len());
        let mut table: Box<[u64]> = vec![0; number_of_entries].into_boxed_slice();

        // querry entire data with for array
        let mut k = 0;
        for i in 0..data.len() {
            for j in i..data.len() {
                table[k] = data[i..=j].iter().copied().min().unwrap();
                k += 1;
            }
        }
        // safe values to array

        Self {
            table: table,
            array_size: data.len(),
        }
    }


    fn space(&self) -> usize {
        std::mem::size_of_val(&*self.table)
    }
        
    fn query(&self, l: usize, r: usize) -> u64 {
        let index = gauss_summation_interval(self.array_size - l, self.array_size) + r - l;
        self.table[index]
    }
}

struct SparseArray {
    sparse_table: Vec<u64>,
    n: usize,
    k: usize,
}

impl<'a> Rmq<'a> for SparseArray {
    fn name() -> String {
        "Sparse Array".to_string()
    }

    fn build(data: &'a [u64]) -> Self {
        let n = data.len();
        
        let k = n.ilog2() as usize;

        let mut sparse_table: Vec<u64> = Vec::with_capacity((n + 1) * (k + 1) + 1 - n);

        for l in 0..=k {
            let interval_size = 0b1 << l;
            for i in 0..n + 1 - interval_size {
                sparse_table.push(data[i..i + interval_size].iter().copied().min().unwrap());

            }
        }


        Self {
            sparse_table,
            n,
            k,
        }
    }

    fn space(&self) -> usize {
        std::mem::size_of_val(&self.sparse_table)
    }
    
    fn query(&self, l: usize, r: usize) -> u64 {
        assert!(l <= r); // requirement
        if l == r { // TODO: check handling and condition in other cases
            return self.sparse_table[l];
        }

        let interval_index = (r - l).ilog2() as usize;

        let mut i: usize = (self.n + 1) * interval_index;
        i -= !(usize::MAX << interval_index);

        let l_true = self.sparse_table[i + l];
        let r_true = self.sparse_table[i + r - (0b1 << interval_index)];
        std::cmp::min(l_true, r_true)
    }
}


struct SegmentTree {
    segment_tree: Vec<Vec<u64>>,
    n: usize,
    k: usize,

}

impl SegmentTree {
        fn recursive_query(&self, l:usize, r:usize, k:usize) -> u64{
        if l == r { // break condition
            return self.segment_tree[k][l];

        }

        let mut next_l = l;
        let mut next_r = r;
        let mut depth_minimum = u64::MAX;
        if l & 0b1 != 0 {
            depth_minimum = self.segment_tree[k][l];
            next_l += 1;
        }
        if r & 0b1 != 0 {
            depth_minimum = std::cmp::min(depth_minimum, self.segment_tree[k][r- 1]);
            next_r -= 1;
        }

        let relevant_diff = std::cmp::min(next_l.trailing_zeros(), next_r.trailing_zeros()) as usize;
        let next_k = k + relevant_diff;
        next_l >>= relevant_diff;
        next_r >>= relevant_diff;
        
        std::cmp::min(depth_minimum, self.recursive_query(next_l, next_r, next_k))
    }

}

impl<'a> Rmq<'a> for SegmentTree {
    fn name() -> String {
        "Segment tree".to_string()
    }

    fn build(data: &'a [u64]) -> Self {
        let n = data.len();
        
        let k = n.ilog2() as usize;

        let mut segment_tree: Vec<Vec<u64>> = vec![Vec::new(); k + 1];

        for l in 0..=k {
            let interval_size = 0b1 << l;
            for i in 0..n / interval_size {
                let block_index = i * interval_size;
                segment_tree[l].push(data[block_index..block_index + interval_size].iter().copied().min().unwrap());

            }
        }


        Self {
            segment_tree,
            n,
            k,
        }
    }

    fn space(&self) -> usize {
        std::mem::size_of_val(self)
    }

    // recursive implementation possible
    fn query(&self, l: usize, r: usize) -> u64 {
        assert!(l <= r); // requirement
        
        self.recursive_query(l, r, 0)
    }
}

struct Blocks {
    segments: SparseArray,
    preffix_table: Vec<u64>,
    suffix_table: Vec<u64>,
    n: usize,
    s: usize,
    block_count: usize,
}

impl<'a> Rmq<'a> for Blocks{
    fn name() -> String {
        "Blocks".to_string()
    }

    fn build(data: &'a [u64]) -> Self {
        let n = data.len();
        let block_size = n.ilog2() as usize;
        let block_count = n / block_size;
        let mut block_minima: Vec<u64> = Vec::with_capacity(block_count);

        // TODO: check index in end and beginning
        for i in 0..block_count {
            block_minima.push(data[i*block_size..(i + 1)*clock_size]);
        }

        let mut prefix_table: Vec<u64> = Vec::with_capacity(n);
        let mut suffix_table: Vec<u64> = Vec::with_capacity(n);

        let mut minimum = u64::MAX;
        for i in 0..n {
            if i % block_size == 0 {
                minimum = u64::MAX;
            }
            suffix_table.push(std::cmp::min(data[i], minimum));
        }

        for i in (0..n).rev() {
            if i % block_size == 0 {
                minimum = u64::MAX;
            }
            prefix_table.push(std::cmp::min(data[i], minimum));
        }

        Self {
            segments: SparseArray::build(block_minima),
            prefix_table: prefix_table,
            suffix_table: suffix_table,
            n: n,
            s: block_size,
            block_count: block_count,
        }
    }

    /// Space usage in bytes.
    fn space(&self) -> usize {
        std::mem::size_of_val(self)
    }

    fn query(&self, l: usize, r: usize) -> u64 {
        let l_block = l / self.s + if l % self.s != 0 { 1 } else { 0 };
        let r_block = r / self.s;

        let pre_minimum = std::cmp::min(self.prefix_table[l], self.segments.query(l_block, r_block));

        std::cmp::min(pre_minimum, self.suffix_table[r]))        
    }
        
   
}

/// The input data.
struct Input {
    data: Vec<u64>,
    queries: Vec<(usize, usize)>,
}

/// Read the given input file.
fn read_input(file: &Path) -> Input {
    let mut input = String::new();
    std::fs::File::open(file)
        .expect("Open input file")
        .read_to_string(&mut input)
        .expect("Read input file");
    let mut vals = input
        .split_ascii_whitespace()
        .map(|s| s.parse::<u64>().unwrap());
    // First line has "{n} {q}"
    let n: usize = vals.next().unwrap() as usize;
    let q: usize = vals.next().unwrap() as usize;
    // Then n lines "{ai}"
    let data = vals.by_ref().take(n).collect();
    // Then q lines "{l} {r}"
    let queries = (0..q)
        .map(|_| (vals.next().unwrap() as usize, vals.next().unwrap() as usize))
        .collect();
    Input { data, queries }
}

/// Bench the given RMQ implementation on the given input, and print the results in CSV format.
fn bench<'a, RMQ: Rmq<'a>>(input: &'a Input) {
    eprint!("{:>10}\t{:>30}\t", input.data.len(), RMQ::name());
    if input.data.len() > RMQ::max_n() {
        eprintln!("skipped");
        return;
    }

    let rmq = RMQ::build(&input.data);
    eprint!("{:>10}\t", rmq.space());
    let start = std::time::Instant::now();
    let mut sum = 0;
    for &(l, r) in &input.queries {
        sum += rmq.query(l, r);
    }
    let elapsed = start.elapsed().as_nanos() as f64 / input.queries.len() as f64;
    println!(
        "{},{},\"{}\",{},{},{}",
        input.data.len(),
        input.queries.len(),
        RMQ::name(),
        rmq.space(),
        sum,
        elapsed
    );
    eprintln!("{:>3}\t{:>8.2}ns/q", sum % 1000, elapsed);
}

fn main() {
    println!("n,q,name,space,sum,time");

    let file_or_dir = PathBuf::from(std::env::args().nth(1).expect("Usage: bench <input_dir>"));

    eprintln!("Reading input from \"{}\" ..", file_or_dir.display());
    let mut inputs = vec![];
    if file_or_dir.is_file() {
        inputs.push(read_input(&file_or_dir));
    } else {
        for entry in file_or_dir.read_dir().expect("Read input directory") {
            if let Ok(file) = entry {
                if Some("in") == file.path().extension().and_then(|s| s.to_str()) {
                    let input = read_input(&file.path());
                    inputs.push(input);
                }
            }
        }
        inputs.sort_by_key(|input| input.data.len());
    }
    for input in inputs {
        bench::<Naive>(&input);
        bench::<LookupTable>(&input);
        bench::<SegmentTree>(&input);
        // TODO: Add other implementations here.
    }
}
