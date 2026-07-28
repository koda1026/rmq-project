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

        // (n + 1 - 2^l) for l in 0..=k
        // (n + 1)(k + 1) - sum(2^l, 0 to =k)
        // (n + 1)(k + 1) - (n - 1)
        // nk + n + k + 1 + 1 - n
        let mut sparse_table: Vec<u64> = Vec::with_capacity((n + 1) * (k + 1) + 1 - n);

        for l in 0..=k {
            let interval_size = 0b1 << l;
            // eprintln!("\nThe current l: {}", l);
            for i in 0..n + 1 - interval_size {
                sparse_table.push(data[i..i + interval_size].iter().copied().min().unwrap());

            }
        }

        // eprintln!("The final output: {}", sparse_table[sparse_table.len() - 1]);

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

        // (n + 1) * l - sum(2^l, 0 to k) (sum of greater and greater powers of 2 is just lsb l
        // bits set to one)
        let mut i: usize = (self.n + 1) * interval_index;
        // eprintln!("The query start index i: {:b} with the interval index: {}, r:{}, l:{}", i, interval_index, r,l);
        i -= !(usize::MAX << interval_index);
        // eprintln!("The query start index i: {:b}", i);
        // (n + 1) * (k + 1) + 1 - n
        // subtract powers of 2 ascending: fill up integer with ones starting on the right side
        // usize::MAX << interval_index and flip

        let l_true = self.sparse_table[i + l];
        let r_true = self.sparse_table[i + r - (0b1 << interval_index)];
        std::cmp::min(l_true, r_true)
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
        bench::<SparseArray>(&input);
        // TODO: Add other implementations here.
    }
}
