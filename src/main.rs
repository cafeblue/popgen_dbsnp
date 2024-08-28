// use clap::{Arg, Command, ArgAction};
use clap::{Args, Parser, Subcommand};

use std::fs::File;
use std::path::Path;
use std::collections::HashMap;
use serde_json::Value;
use std::io::{self, BufRead, BufReader};
use bzip2::read::BzDecoder;
use itertools::Itertools;

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn fac(n: u64) -> u64 {
    if n < 2 {
        return 1;
    }
    else {
        return n * fac(n-1);
    }
}

fn binom(n: u64, k: u64) -> f64 {
    if n < k {
        return 0.0;
    }
    else {
        return fac(n) as f64 / (fac(k) as f64 * fac(n-k) as f64);
    }
}

fn fn_a_1(all_nb: u64) -> f64 {
    let mut tmp_a_1: f64 = 0.0;
        for m in 1..(all_nb+1) {
            tmp_a_1 += 1.0 / m as f64;
        }
        return tmp_a_1;
}

fn tajd_cstes(nb: u64) -> (f64, f64) {
	let mut a1: f64 = 0.0;
	let mut a2: f64 = 0.0;
	for n in 1..nb {
		a1 += 1.0 / n as f64;
		a2 += 1.0 / (n.pow(2) as f64);
	}

	let b1: f64 = (nb + 1) as f64 / ((3 * (nb - 1)) as f64);
	let b2: f64 = ((2 * (nb.pow(2) + nb + 3)) as f64 )/ ((9 * nb * (nb - 1)) as f64);
	let c1: f64 = b1 - (1.0 / a1);
	let c2: f64 = b2 - ((nb + 2) as f64 / (a1 * nb as f64)) + (a2 / a1.powi(2));
	let e1: f64 = c1 / a1;
	let e2: f64 = c2 / (a1.powi(2) + a2);
	return (e1, e2);
}

fn bai(sh: u64) -> f64 {
    let mut beh: f64 = 0.0;
    for gna in 1..sh {
        beh += 1.0/(gna.pow(2) as f64);
    }
    return beh;
}

 
fn hnorm(enn: u64, ess: u64, teta: f64) -> f64 {
    let ft: f64 = teta * ( enn as f64 - 2.0) / (6.0 * (enn as f64 - 1.0));
    let mut aah: f64 = 0.0;
    for gnii in 1..enn {
        aah += 1.0 / gnii as f64;
    }

    let tetacarre: f64 = ess as f64 * (ess as f64 - 1.0) / ( aah.powi(2) as f64 + bai(enn) );
    let thrd: f64 = 18.0 * enn.pow(2) as f64 * ( 3.0 * enn as f64 + 2.0) * bai(enn+1) - (88.0 * enn.pow(3) as f64 + 9.0 * enn.pow(2) as f64 - 13.0 * enn as f64 + 6.0 ) / ( 9.0 * enn as f64 * (enn-1).pow(2) as f64);
    let vari: f64 = thrd * tetacarre + ft;
    let mut sqv: f64 = 0.0;
    if vari > 0.0 {
        sqv = vari.sqrt();
    }
    return sqv;
}

fn etest_under(enn: u64, ess: u64, teta: f64) -> f64 {
    let mut beh: f64 = 0.0;
    for gna in 1..enn {
        beh += 1.0/(gna.pow(2) as f64);
    }

    let mut aah: f64 = 0.0;
    for gnii in 1..enn {
        aah += 1.0 / gnii as f64;
    }


    let ft: f64 = teta * ( enn as f64 / (enn as f64 - 1.0) / 2.0 - 1.0 / aah );

    let tetacarre: f64 = ess as f64 * (ess as f64 - 1.0) / ( aah.powi(2) as f64 + bai(enn) );
    let thrd: f64 = beh / aah.powi(2) + (enn as f64 / (enn as f64 - 1.0)).powi(2) * beh * 2.0 - 2.0 * (enn as f64 * beh - enn as f64 + 1.0) / (enn as f64 - 1.0) / aah - (enn as f64 * 3.0 + 1.0) / (enn as f64 - 1.0); 
    let vari: f64 = thrd * tetacarre + ft;
    let mut sqv: f64 = 0.0;
    if vari > 0.0 {
        sqv = vari.sqrt();
    }
    return sqv;
}

fn dbsnp2af(dbsnp_file: String, gtf_file: String, study_name: String, bin_size: u16, exon_only: bool) {
  let mut exon_posset: HashMap<String, (HashMap<u32, u16>, u16)> = HashMap::new();
  eprintln!("Study_Name: {}", study_name);
  eprintln!("bin size: {}", bin_size);
  eprintln!("Use exon only: {}", exon_only);

  // Read in values
  if let Ok(lines) = read_lines(gtf_file) {
    for line in lines {
      if let Ok(ip) = line {
        let spl: Vec<&str> = ip.split("\t").collect();
        if spl[2] != "exon" && exon_only{
          continue;
        }
        if spl[2] != "gene" && !exon_only{
          continue;
        }
        let startp = spl[3].parse::<u32>().unwrap();
        let stopp = spl[4].parse::<u32>().unwrap();
        let mut currexon: HashMap<u32, u16> = HashMap::new();
        for i in startp..(stopp+1) {
          currexon.insert(i, 0);
        }
        let desc = spl[8];
        let desc_parts: Vec<&str> = desc.split("; ").collect();
        for part in desc_parts {
          let part_parts: Vec<&str> = part.split(" ").collect();
          if part_parts[0] == "gene_name" {
            let gene_id = part_parts[1].replace("\"", "");
            if !exon_posset.contains_key(&gene_id) {
              let tup: (HashMap<u32, u16>, u16) = (currexon.clone(), 0);
              exon_posset.insert(gene_id, tup.clone());
            }
            else {
              let (tmp, tmpbin) = exon_posset.get(&gene_id).unwrap();
              currexon.extend(tmp.clone().into_iter());
              let tup: (HashMap<u32, u16>, u16) = (currexon.clone(), tmpbin.clone());
              exon_posset.insert(gene_id, tup.clone());
            }
          }
        }
      }
    }
  }
  
  if bin_size > 0 {
    let mut exon_posset_sorted: HashMap<String, (HashMap<u32, u16>, u16)> = HashMap::new();
    for gene in exon_posset.keys() {
      let mut cbin: u16 = 0 ;
      let mut cpos: u16 = 0 ;
      let mut currexon: HashMap<u32, u16> = HashMap::new();
      let (tmpposset, _rmme) = exon_posset.get(gene).unwrap(); 
      for positions in tmpposset.keys().sorted() {
        if cpos < &bin_size * (cbin + 1) {
          currexon.insert(positions.clone(), cbin.clone());
        }
        else {
          cbin += 1;
          currexon.insert(positions.clone(), cbin.clone());
        }
        cpos += 1;
      }
      let tup: (HashMap<u32, u16>, u16) = (currexon.clone(), cbin.clone()); 
      exon_posset_sorted.insert(gene.to_string(), tup.clone());
    }
    exon_posset = exon_posset_sorted;
  }
  
  let mut gene_tmppi: HashMap<String, HashMap<PiTmp, u16>> = HashMap::new();
  let file = File::open(dbsnp_file).unwrap();
  let reader = BufReader::new(BzDecoder::new(file));
  for line in reader.lines() {
    let l = line.unwrap();
    if !l.contains( &study_name) {
      continue;
    }
    let v: Value = serde_json::from_str(&l).unwrap();

    let mut a_count: u64 = 0;
    let mut c_count: u64 = 0;
    let mut freq_idx: usize = 127;
    let mut tmp_pos: u32 = 0;
    let mut on_exon: bool = false;

    'marker: for idx in 0..v["primary_snapshot_data"]["placements_with_allele"].as_array().unwrap().len() {
      if v["primary_snapshot_data"]["placements_with_allele"][idx]["is_ptlp"].as_bool().unwrap() {
        for idx1 in 0..v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"].as_array().unwrap().len() {
          if v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["deleted_sequence"].as_str().unwrap() == v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["inserted_sequence"].as_str().unwrap() {
            tmp_pos = v["primary_snapshot_data"]["placements_with_allele"][idx]["alleles"][idx1]["allele"]["spdi"]["position"].as_u64().unwrap() as u32;
            freq_idx = idx1;
            break 'marker;
          }
        }
      }
    }

    let gene_symbol_option = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["locus"].as_str();

    let mut bintmp = 0;
    let mut tbintmp = 0;
    if let Some(gene_symbol) = gene_symbol_option {
    // gene_symbol is not None, you can use it here
      if exon_posset.contains_key(&gene_symbol.to_string()) {
        let (tmp , tmpbin)= exon_posset.get(&gene_symbol.to_string()).unwrap();
        if tmp.contains_key(&tmp_pos) {
          on_exon = true;
          bintmp = tmp.get(&tmp_pos).unwrap().clone();
          tbintmp = tmpbin.clone();
        }
        else {
          continue;
        }
      } 
    } else {
      continue;
    }        

    if freq_idx == 127 {
      eprintln!("failed to get freq_idx: {}", v["refsnp_id"]);
      continue;
    }
    let gene_symbol = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["locus"].as_str().unwrap();
    // let vt1 = v["primary_snapshot_data"]["allele_annotations"][freq_idx]["assembly_annotation"][0]["genes"][0]["rnas"][0]["sequence_ontology"][0]["name"].to_string();
    for idx in 0..v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"].as_array().unwrap().len() {
      //if args.exac {
      if v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["study_name"].as_str().unwrap() == study_name{
        a_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["allele_count"].as_u64().unwrap() as u64;
        c_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["total_count"].as_u64().unwrap() as u64;
      }
      //}
      //else {
      //    a_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["allele_count"].as_u64().unwrap() as u64;
      //    c_count += v["primary_snapshot_data"]["allele_annotations"][freq_idx]["frequency"][idx]["total_count"].as_u64().unwrap() as u64;
      //}
    }
    // let vt2: String = v["primary_snapshot_data"]["variant_type"].to_string();
    if a_count == 100 || c_count < 100 || a_count == c_count || !on_exon {
      continue;
    }
    let tmp_pos: PiTmp = PiTmp {total: c_count, refcount: a_count, bin: bintmp, totalbin: tbintmp};

    if let Some(tmp) = gene_tmppi.get_mut(&gene_symbol.to_string()) {
      if tmp.contains_key(&tmp_pos) {
        let tmp2 = tmp.get(&tmp_pos).unwrap() + 1;
        tmp.insert(tmp_pos, tmp2);
      }
      else {
        tmp.insert(tmp_pos, 1);
      }
    }
    else {
      gene_tmppi.insert(gene_symbol.to_string(), HashMap::from([(tmp_pos, 1),]));
    }
  };

  for gene in gene_tmppi.keys() {
    let tmp = gene_tmppi.get(gene).unwrap();
    for pos in tmp.keys() {
      println!("{}\t{}\t{}\t{}\t{}\t{}", gene, pos.refcount, pos.total, tmp.get(pos).unwrap(), pos.bin, pos.totalbin);
    }
  }
}

fn popgen(af_file: String, sfs_file: String) {
  println!("af_file: {}", af_file);
  /*let mut sfs: HashMap<u32, u32> = HashMap::new();
  let mut gene_tmppi: HashMap<String, HashMap<PiTmp, u16>> = HashMap::new();
  let file = File::open(sfs_file).unwrap();
  let reader = BufReader::new(file);
  for line in reader.lines() {
    let l = line.unwrap();
    let spl: Vec<&str> = l.split("\t").collect();
    sfs.insert(spl[0].parse::<u32>().unwrap(), spl[1].parse::<u32>().unwrap());
  }
  let file = File::open(af_file).unwrap();
  let reader = BufReader::new(file);
  for line in reader.lines() {
    let l = line.unwrap();
    let spl: Vec<&str> = l.split("\t").collect();
    let tmp_pos: PiTmp = PiTmp {total: spl[2].parse::<u32>().unwrap(), refcount: spl[1].parse::<u32>().unwrap() , bin: spl[4].parse::<u32>().unwrap(), totalbin: spl[5].parse::<u32>().unwrap()};
    gene_tmppi.insert(spl[0].parse::<String>().unwrap(), HashMap::from([(tmp_pos, spl[3].parse::<u16>().unwrap()),]));
  }
  
  for gene in gene_tmppi.keys() {
    let tmp = gene_tmppi.get(gene).unwrap();

    let (poslist, binnumber) = exon_posset.get(gene).unwrap();
    let totaldet = poslist.len() as u32;
    let mut det: Vec<f64> = vec![0.0; binnumber.clone() as usize + 1];
    if args.bin > 0 {
      for i in 0..(binnumber.clone() as usize) {
        det[i] = args.bin.clone() as f64;
      }
      det[binnumber.clone() as usize] = (totaldet % args.bin) as f64;
    }
    else {
        det[0] = totaldet as f64;
    }
    let mut all_nb: Vec<u64> = vec![1; binnumber.clone() as usize + 1];
    let mut pi_tmp: Vec<u64> = vec![0; binnumber.clone() as usize + 1];
    let mut pi_min_tmp: Vec<u64> = vec![0; binnumber.clone() as usize + 1];
    let mut h_tmp: Vec<u64> = vec![0; binnumber.clone() as usize + 1];
    let mut hn_tmp: Vec<u64> = vec![0; binnumber.clone() as usize + 1];
    let mut seg: Vec<u64> = vec![0; binnumber.clone() as usize + 1];
    for pos in tmp.keys() {
        let tmpbinid = pos.bin as usize;
        if pos.total > all_nb[tmpbinid] {
            all_nb[tmpbinid] = pos.total.clone();
        }
        pi_tmp[tmpbinid] += pos.refcount * (pos.total - pos.refcount) * (tmp.get(pos).unwrap().clone() as u64);
        pi_min_tmp[tmpbinid] += 1 * (pos.total - 1) * (tmp.get(pos).unwrap().clone() as u64);
        h_tmp[tmpbinid] +=  pos.refcount.pow(2) * (tmp.get(pos).unwrap().clone() as u64);
        hn_tmp[tmpbinid] += pos.refcount * (tmp.get(pos).unwrap().clone() as u64) ;
        seg[tmpbinid] += tmp.get(pos).unwrap().clone() as u64;
        //println!("{}\t{}\t{}\t{}\t{}", gene, pos.refcount, pos.total, tmp.get(pos).unwrap(), tmpbinid);
    }
	//eprintln!("all_nb: {}", all_nb);
    // calculate a_1
    for i in 0..(binnumber.clone() as usize + 1) {
        let mut tajd: f64 = 0.0; 
        let mut tajd_min: f64 = 0.0; 
        let a_1 = fn_a_1(all_nb[i]);

        let pi = pi_tmp[i] as f64 / ((all_nb[i] as f64 * (all_nb[i] as f64 - 1.0)) / 2.0) / det[i] ;
        let pi_min = pi_min_tmp[i] as f64 / ((all_nb[i] as f64 * (all_nb[i] as f64 - 1.0)) / 2.0) / det[i] ;
        let theta = seg[i] as f64 / a_1 / det[i];

        let mut d: f64 = 0.0;
        let mut dmin: f64 = 0.0;
        let (e1, e2) = tajd_cstes(all_nb[i]);
        let p1p2: f64 = e1 * seg[i] as f64 + e2 * seg[i] as f64 * (seg[i] as f64 - 1.0);
        if pi > 0.0 && theta > 0.0 {
          	d = pi - theta;
      	    dmin = pi_min - theta;
            dmin = dmin.abs();
        }

        if p1p2 > 0.0 {
            tajd = d * det[i] / p1p2.sqrt();
            tajd_min = dmin * det[i] / p1p2.sqrt();
        }

        let hache = h_tmp[i] as f64/ ((all_nb[i] as f64 * (all_nb[i] as f64 - 1.0)) / 2.0);
        let final_h = pi - hache;

        let elle = hn_tmp[i] as f64 / (all_nb[i] as f64 - 1.0);
        let mut final_hnorm: f64 = 0.0;
        let gun = hnorm(all_nb[i], seg[i], theta);
        if gun > 0.0 {
            final_hnorm = (pi - elle) / gun 
        }

        let mut final_e: f64 = 0.0;
        let eun: f64 = etest_under(all_nb[i], seg[i], theta);
        if eun > 0.0 {
            final_e = (elle - theta) / eun;
        }


        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", gene, all_nb[i], totaldet, i, det[i], pi, theta, tajd, tajd/tajd_min, final_h, final_hnorm, final_e);
    }
  }*/
}

/*
for my $spl (2..$all_nb){
	for my $u (1..$all_nb-1){		## foreach existing sample size - so for each potential nb of spls in SFS
		for my $v (0..$spl){	
			$SFSbins{$spl}{$u}{$v}=( ( &binom($all_nb - $u, $spl - $v) * &binom( $u, $v ) ) ) / &binom( $all_nb, $spl ); 
		}
	} */


fn sfs(samplesize: u64) {
  for m in 2..(samplesize+1) {
    for n in 1..samplesize {
      for v in 0..(m+1) {
        let mut tmp: f64 = 0.0;
        tmp = binom(samplesize - n, m - v) * binom(n, v) / binom(samplesize, m);
        println!("{}\t{}\t{}\t{}", m, n, v, tmp);
      }
    }  

  }
}

#[derive(Parser, Debug)]
#[clap(name = "popgen_dbsnp", author = "Wei Wang", version = "0.1.0", about = "generate popgen stat based on dbsnp file.", long_about = None)]
struct App {
    #[clap(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Parse dbsnp json file and gtf file, generate a table of allic frequency.
    Parse2af(ParseDBSNP),

    /// Calculate population genetics statistics based on the allele frequency table.
    Popgen(CalculatePOPGen),

    /// Calculate SFS based on the sample size.
    Sfs(SFSbin),
}

#[derive(Args, Debug)]
struct ParseDBSNP {
    /// dbsnp json file
    #[arg(long)]
    dbsnp: String,

    /// gtf file
    #[arg(long)]
    gtf: String,

    /// study name
    #[arg(long)]
    studyname: String,

    /// bin size
    #[arg(long)]
    bin: u16,

    /// use exon only
    #[arg(long)]
    exononly: bool,
}

#[derive(Args, Debug)]
struct CalculatePOPGen {
    /// sfs file
    #[arg(long)]
    sfstable: String,

    /// af file
    #[arg(long)]
    aftable: String,
}

#[derive(Args, Debug)]
struct SFSbin {
    /// sfs file
    #[arg(long)]
    samplesize: u64,
}


#[derive(Debug, Eq, PartialEq, Hash)]
struct PiTmp {
    total: u64,
    refcount: u64,
    bin: u16,
    totalbin: u16,
}
fn main() {
  let args = App::parse();
  match &args.command {
    Command::Parse2af(name) => {
        dbsnp2af(name.dbsnp.clone(), name.gtf.clone(), name.studyname.clone(), name.bin.clone(), name.exononly.clone());
    }
    Command::Popgen(name) => {
        popgen(name.sfstable.clone(), name.aftable.clone());
    }
    Command::Sfs(name) => {
        sfs(name.samplesize.clone());
    }
  }
}


