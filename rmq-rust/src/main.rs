use std::{
    collections::{HashMap, VecDeque},
    fs,
    hash::Hash,
    io::Read,
    mem::{size_of, size_of_val},
    path::{Path, PathBuf},
    cmp,
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
    sparse_table: Vec<Vec<u64>>,
    n: usize,
}

impl<'a> Rmq<'a> for SparseArray {
    fn name() -> String {
        "Sparse Array".to_string()
    }

    fn max_n() -> usize {
        // NOTE: Do not use this for the improved implementations!
        10_000
    }
    fn build(data: &'a [u64]) -> Self {
        let n = data.len();
        
        let k = n.ilog2() as usize;

        let mut sparse_table: Vec<Vec<u64>> = Vec::with_capacity(k + 1);

        for l in 0..=k {
            let interval_size = 0b1 << l;
            let mut interval_table: Vec<u64> = Vec::with_capacity(n + 1 - interval_size);
            for i in 0..n + 1 - interval_size {
                interval_table.push(data[i..i + interval_size].iter().copied().min().unwrap());

            }

            sparse_table.push(interval_table);
        }


        Self {
            sparse_table,
            n,
        }
    }

    fn space(&self) -> usize {
        // space of the level list
        let mut total_size = size_of_val(&self.sparse_table);
        total_size += size_of::<Vec<u64>>() * self.sparse_table.capacity();

        // space of each individual levels
        for l in 0..self.sparse_table.len() {
            total_size += size_of::<u64>() * self.sparse_table[l].capacity();
        }

        return total_size;
    }
    
    fn query(&self, l: usize, r: usize) -> u64 {
        assert!(l <= r); // requirement
        let depth = ((r + 1) - l).ilog2() as usize;

        let l_min = self.sparse_table[depth][l];
        let r_min = self.sparse_table[depth][(r + 1) - (0b1 << depth)];
        cmp::min(l_min, r_min)
    }
}


struct SegmentTree {
    segment_tree: Vec<Vec<u64>>,
    k: usize,
}

impl<'a> Rmq<'a> for SegmentTree {
    fn name() -> String {
        "Segment tree".to_string()
    }

    fn max_n() -> usize {
        // NOTE: Do not use this for the improved implementations!
        10_000
    }

    fn build(data: &'a [u64]) -> Self {
        let n = data.len();
        
        let k = n.ilog2() as usize;

        let mut segment_tree: Vec<Vec<u64>> = Vec::with_capacity(k + 1);

        for l in 0..=k {
            let interval_size = 0b1 << l;
            let mut interval_table = Vec::with_capacity(n / k);
            for i in 0..n / interval_size {
                let block_index = i * interval_size;
                interval_table.push(data[block_index..block_index + interval_size].iter().copied().min().unwrap());

            }
            segment_tree.push(interval_table);
        }


        Self {
            segment_tree,
            k,
        }
    }

    fn space(&self) -> usize {
        // space of the level list
        let mut total_size = size_of_val(&self.segment_tree);
        total_size += size_of::<Vec<u64>>() * self.segment_tree.capacity();

        // space of each individual levels
        for l in 0..self.segment_tree.len() {
            total_size += size_of::<u64>() * self.segment_tree[l].capacity();
        }

        total_size += size_of_val(&self.k);

        return total_size;
    }
 
    // recursive implementation possible
    fn query(&self, l: usize, r: usize) -> u64 {
        let mut new_l = l;
        let mut new_r = r + 1;

        let mut total_minimum = u64::MAX;
        let mut level = 0;
        while new_l != new_r {
            if new_l & (0b1 << level) != 0{
                let row = &self.segment_tree[level];
                total_minimum = cmp::min(row[new_l >> level], total_minimum);
                new_l += 0b1 << level;
            }
            if new_r & (0b1 << level) != 0 {
                total_minimum = cmp::min(self.segment_tree[level][(new_r >> level) - 1], total_minimum);
                new_r -= 0b1 << level;
            }
            level += 1;
        }
        
        return total_minimum;
    }
}

struct Blocks<'a> {
    blocks: SparseArray,
    prefix_table: Vec<u64>,
    suffix_table: Vec<u64>,
    data: &'a [u64],
    cartesian_tree: CartesianTree,
    s: usize,
    block_count: usize,
    pre_calc: bool,
}

impl<'a> Rmq<'a> for Blocks<'a>{
    fn name() -> String {
        "Blocks".to_string()
    }

    fn max_n() -> usize {
        // NOTE: Do not use this for the improved implementations!
        10_000
    }
    fn build(data: &'a [u64]) -> Self {
        let pre_calc = false;
        let with_cartesian_tree = true;
        let n = data.len();
        let block_size = n.ilog2() as usize;
        let block_count = n / block_size;
        let mut block_minima:Vec<u64> = Vec::with_capacity(block_count);

        // TODO: check index in end and beginning
        for i in 0..block_count {
            block_minima.push(data[i*block_size..(i + 1)*block_size].iter().copied().min().unwrap());
        }

        let mut prefix_table: Vec<u64>;
        let mut suffix_table: Vec<u64>;
        let cartesian_tree: CartesianTree = CartesianTree::build(data, block_size);

        if pre_calc {
            prefix_table = vec![u64::MAX; n];
            suffix_table = Vec::with_capacity(n);
            let mut minimum = u64::MAX;
            for i in 0..n {
                if i % block_size == 0 {
                    minimum = u64::MAX;
                }
                minimum = cmp::min(data[i], minimum);
                suffix_table.push(minimum);
            }

            minimum = u64::MAX;

            for i in (0..n).rev() {
                if i % block_size == 0 {
                    minimum = u64::MAX;
                }
                minimum = cmp::min(data[i], minimum);
                prefix_table[i] = minimum;
            }
        } else {
            prefix_table = Vec::with_capacity(0);
            suffix_table = Vec::with_capacity(0);
        }

        Self {
            blocks: SparseArray::build(&block_minima),
            prefix_table: prefix_table,
            suffix_table: suffix_table,
            data: data,
            cartesian_tree: cartesian_tree,
            s: block_size,
            block_count: block_count,
            pre_calc: pre_calc,
        }
    }

    /// Space usage in bytes.
    fn space(&self) -> usize {
        self.blocks.space() +
        self.prefix_table.capacity() * size_of::<u64>() + size_of_val(&self.prefix_table) +
        self.suffix_table.capacity() * size_of::<u64>() + size_of_val(&self.suffix_table) +
        size_of_val(self.data) +
        size_of_val(&self.s) +
        size_of_val(&self.block_count) +
        size_of_val(&self.pre_calc)
    }

    fn query(&self, l: usize, r: usize) -> u64 {
        let mut l_block = l / self.s;
        if l % self.s != 0 {
            l_block += 1;
        }
        let mut r_block = (r / self.s);

        let mut total_minimum;
        if l_block < r_block {
            assert!(r_block > 0);
            let left_minimum;
            let right_minimum;
            if self.pre_calc {
                left_minimum = self.prefix_table[l];
                right_minimum = self.suffix_table[r];
            } else {
                let old_left_minimum = self.data[l..=l_block * self.s].iter().copied().min().unwrap();
                let old_right_minimum = self.data[r_block * self.s..=r].iter().copied().min().unwrap();
                left_minimum = self.data[self.cartesian_tree.prefix_query(l)];
                if r >= self.data.len() / self.s * self.s {
                    right_minimum =  self.data[r_block * self.s..=r].iter().copied().min().unwrap();
                } else {
                    right_minimum = self.data[self.cartesian_tree.suffix_query(r)];
                }
            }
            r_block -= 1; // inclusive index
            let block_minimum = self.blocks.query(l_block, r_block);
            total_minimum = cmp::min(cmp::min(left_minimum, right_minimum), block_minimum);
        } else { // case where l and r are in the same block
            total_minimum = self.data[l..=r].iter().copied().min().unwrap();
        }

        return total_minimum;
    }
}

struct TreeNumber{
    number: VecDeque<usize>,
    bit_count: usize,
}

impl TreeNumber{
    fn add_child_value(&mut self, child_number: &Self) {
        // 1. appropriately shift value
        let size_of_block = std::mem::size_of_val(&usize::MAX) * 8; // in bits not in bytes
        let additional_blocks = child_number.bit_count / size_of_block;
        let left_over_bits = child_number.bit_count % size_of_block;
        // 1.1 move bits in vector(copy elements, then shift)
        if (self.bit_count % size_of_block + left_over_bits) > size_of_block {
            self.number.insert(self.number.len(), 0); // push to end
        }

        let mut buffer = 0;
        let mut prev_bits = 0;
        for i in 0..self.number.len(){
            buffer = self.number[i];
            self.number[i] <<= left_over_bits;
            self.number[i] |= prev_bits;
            prev_bits = buffer >> size_of_block - left_over_bits;
        }

        // 1.2 adjust size of vector
        for i in 0..additional_blocks {
            self.number.insert(0, 0);
        }

        // 2. add child number
        for i in 0..child_number.number.len() {
            self.number[i] |= child_number.number[i];
        }

        // 3. update meta data
        self.bit_count += child_number.bit_count;
    }

    fn add_dead_end(&mut self) {
        self.add_child_value(&Self { number: VecDeque::from([0b0]), bit_count: 1 });
    }

    fn get_number(&self) -> VecDeque<usize> {
        self.number.clone()
    }

    // preorder bit encoding with leaves encoded as zeros
    fn build(block: &[u64]) -> TreeNumber {
        let mut result = TreeNumber{ number: VecDeque::from([0b1]), bit_count:1 };

        let mut min_index = 0;
        let mut block_min = u64::MAX;
        for i in 0..block.len() {
            if block[i] < block_min {
                min_index = i;
                block_min = block[i];
            }
        }


        if min_index > 0 {
            let left_child = TreeNumber::build(&block[0..min_index]);
            result.add_child_value(&left_child);
        } else {
            result.add_dead_end(); // Add zero to end of tree representation (no child)
        }

        if min_index < block.len() - 1 {
            let right_child = TreeNumber::build(&block[min_index + 1..block.len()]);
            result.add_child_value(&right_child);
        } else {
            result.add_dead_end(); // Add zero to end of tree representation (no child)
        }

        return result;
    }
}

struct CartesianTree {
    tree_numbers: Vec<VecDeque<usize>>,
    tree_lookup_tables: HashMap<VecDeque<usize>,(Vec<usize>, Vec<usize>)>,
    block_size: usize,
}

impl CartesianTree {
    fn build_lookup_table_old(block: &[u64]) -> Vec<Vec<usize>> {
        let mut lookup_table:Vec<Vec<usize>> = Vec::with_capacity(block.len());

        for l in 0..block.len() {
            let mut row = Vec::with_capacity(block.len() - l);
            for r in l..block.len() {
                row.push(block[l..=r].iter().copied().position(|x| block[l..=r].iter().all(|&y| x <= y)).unwrap() + l);
            }
            lookup_table.push(row);
        }

        return lookup_table;
    }

    fn build_lookup_table(block: &[u64]) -> (Vec<usize>, Vec<usize>) {
        let mut prefix_table = vec![usize::MAX; block.len()];
        let mut suffix_table = Vec::with_capacity(block.len());

        let mut minimum = u64::MAX;
        let mut min_index = 0;
        for i in 0..block.len() {
            if block[i] < minimum {
                min_index = i;
                minimum = block[i];
            }
            suffix_table.push(min_index);
        }

        minimum = u64::MAX;
        min_index = 0;

        for i in (0..block.len()).rev() {
            if block[i] < minimum {
                min_index = i;
                minimum = block[i];
            }
            prefix_table[i] = min_index;
        }

        (prefix_table, suffix_table)
    }
   
    fn build(data: &[u64], block_size: usize) -> Self {
        // let block_size = (data.len().ilog2() / 4) as usize;
        let block_count = data.len() / block_size;
        let mut tree_numbers:Vec<VecDeque<usize>> = Vec::with_capacity(block_count);
        let mut tree_lookup_tables:HashMap<VecDeque<usize>, (Vec<usize>, Vec<usize>)> = HashMap::new();
        
        for block_index in 0..block_count {
            let block = &data[block_index * block_size..(block_index + 1) * block_size];
            let tree_number = TreeNumber::build(&block);
            tree_numbers.push(tree_number.get_number());

            if !tree_lookup_tables.contains_key(&tree_number.get_number()) {
                let lookup_table = CartesianTree::build_lookup_table(&block);
                tree_lookup_tables.insert(tree_number.get_number(), lookup_table);
            }
        }

        Self {
            tree_numbers,
            tree_lookup_tables,
            block_size,
        }
    }

    fn prefix_query(&self, l: usize) -> usize {
        let mut block_index = l / self.block_size;

        let tree_number = &self.tree_numbers[block_index];
        let lookup_table = &self.tree_lookup_tables[&tree_number].0;
        let block_minimum_index = lookup_table[l % self.block_size];

        block_minimum_index + block_index * self.block_size
    }

    fn suffix_query(&self, r: usize) -> usize {
        let block_index = r / self.block_size;

        let tree_number = &self.tree_numbers[block_index];
        let lookup_table = &self.tree_lookup_tables[&tree_number].1;
        let block_minimum_index = lookup_table[r % self.block_size];

        block_minimum_index + block_index * self.block_size
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
    let mut sum:u64 = 0;
    for &(l, r) in &input.queries {
        let minimum = rmq.query(l, r);
        sum = sum.wrapping_add(minimum);
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
    //println!("n,q,name,space,sum,time");

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
    let mut query_monitor:HashMap<(usize,usize), u64> = HashMap::new();
    for input in inputs {
        bench::<Naive>(&input);
        //bench::<LookupTable>(&input);
        // bench::<SparseArray>(&input);
        // bench::<SegmentTree>(&input);
        bench::<Blocks>(&input);
        // bench::<CartesianTree>(&input);
        // TODO: Add other implementations here.
    }
}
