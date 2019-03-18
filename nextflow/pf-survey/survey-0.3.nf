log.info """
=======================================================
                     PathFinder
                Survey Pipeline v0.3
               (QC, Assembly, Typing)
=======================================================

run_id          =    $params.run_id
out_dir         =    $params.out_dir
file_dir        =    $params.file_dir
file_pattern    =    $params.file_pattern
file_extension  =    $params.file_extension
resource_dir    =    $params.resource_dir

leading         =    $params.leading
trailing        =    $params.trailing
min_len         =    $params.min_len
sliding_window  =    $params.sliding_window
adapters        =    $params.adapters
kraken_db       =    $params.minikraken
assembler       =    $params.assembler
depth           =    $params.depth
reference       =    $params.reference

mash            =    $params.mash
mykrobe         =    $params.mykrobe
kleborate       =    $params.kleborate

=======================================================
=======================================================
"""


fastq = Channel
		.fromFilePairs("${params.file_dir}/${params.file_pattern}${params.file_extension}", flat: true)
		.ifEmpty { exit 1, "Input read files could not be found." }

if (params.genus) {
    use_genus = "--usegenus"
} else {
    use_genus = ""
}

process Trimmomatic {

    label "survey"
	label "medium"

	tag { id }

	input:
	set id, file(forward), file(reverse) from fastq

	output:
	set id, file("${id}_1.fq.gz"), file("${id}_2.fq.gz") into kraken, shovill, mykrobe

	"""
	trimmomatic PE ${forward} ${reverse} -threads $task.cpus $params.phred -baseout ${id}.fq.gz \
	ILLUMINACLIP:${params.resource_dir}/${params.adapters}:2:30:10:3:TRUE LEADING:${params.leading} TRAILING:${params.trailing} \
	SLIDINGWINDOW:${params.sliding_window} MINLEN:${params.min_len}

	mv ${id}_1P.fq.gz ${id}_1.fq.gz
	mv ${id}_2P.fq.gz ${id}_2.fq.gz
	"""

}

process Shovill {

  label "survey"
	label "high"

	tag { id }
	publishDir "${params.out_dir}/skesa", mode: "copy"

	input:
	set id, file(forward), file(reverse) from shovill.filter{ it[1].size()>10000 & it[2].size()>10000 } // Failed / empty file can be up to 10kb

	output:
	set id, file("${id}.fasta") into mash, abricate, mlst, prokka, kleborate

	"""
	shovill --R1 ${forward} --R2 ${reverse} --cpus $task.cpus --ram $task.memory \
  --depth $params.depth --assembler $params.assembler --outdir $id --force
  mv ${id}/contigs.fa ${id}.fasta
	"""

}


process Kraken2 {

    label "survey"
  	label "high"

    tag { id }
    publishDir "${params.out_dir}/kraken", mode: "copy"

    input:
    set id, file(forward), file(reverse) from kraken.filter{ it[1].size()>10000 & it[2].size()>10000 }

    output:
    set id, file("${id}.report")

    """

    kraken2 --db ${params.resource_dir}/${params.minikraken} --threads $task.cpus --output "-" \
    --fastq-input --gzip-compressed --paired --report ${id}.report --use-names ${forward} ${reverse}
    """

}

if (params.mykrobe) {
    process MykrobePredictor {

        label "mykrobe"
          label "medium"

        tag { id }
        publishDir "${params.out_dir}/mykrobe", mode: "copy"

        input:
        set id, file(forward), file(reverse) from mykrobe.filter{ it[1].size()>10000 & it[2].size()>10000 }

        output:
        set id, file("${id}.json")

        """
        mykrobe predict --threads $task.cpus ${id} $params.mykrobe --seq $forward $reverse > ${id}.json
        """

    }
}

process Abricate {

  label "survey"
	label "minimum"

	// Assembly-based typing with T. Seemann's Abricate

	tag { id }

	publishDir "${params.out_dir}/abricate", mode: "copy"

	input:
	set id, file(assembly) from abricate.filter{ it[1].size()>10000 }
	each abricate_db from params.abricate_dbs

	output:
	file("${abricate_db}/${id}.tab")

	"""
	mkdir -p ${abricate_db}
	abricate --db ${abricate_db} ${assembly} > ${abricate_db}/${id}.tab
	"""

}

process MLST {

  label "survey"
	label "minimum"

	// Assembly-based MLST

	tag { id }
	publishDir "${params.out_dir}/mlst", mode: "copy"

	input:
	set id, file(assembly) from mlst.filter{ it[1].size()>10000 }

	output:
	file("${id}.tab")

	"""
	mlst $assembly > ${id}.tab
	"""

}

if (params.mash){
    process Mash {

        label "survey"
        label "minimum"

        // MinHash distance against global reference:

        tag { id }
        publishDir "${params.out_dir}/mash", mode: "copy"

        input:
        set id, file(assembly) from mash.filter{ it[1].size()>10000 }

        output:
        file("${id}.mash.tab")

        """
        mash dist $assembly $params.resource_dir/$params.reference > ${id}.mash.tab
        """
    }
}

process Prokka {

  label "survey"
  label "medium"

  // MinHash distance against global reference:

  tag { id }
  publishDir "${params.out_dir}/prokka", mode: "copy"

  input:
  set id, file(assembly) from prokka.filter{ it[1].size()>10000 }

  output:
  file("${id}/${id}.gff")

  """
  prokka --outdir ${id} --force --prefix ${id} --kingdom $params.kingdom --genus $params.genus \
  $use_genus --fast --cpus $task.cpus $assembly
  """

}
if (params.kleborate) {

  process Kleborate {

    label "kleb"
    label "minimal"

    // MinHash distance against global reference:

    tag { id }
    publishDir "${params.out_dir}/kleborate", mode: "copy"

    input:
    set id, file(assembly) from kleborate.filter{ it[1].size()>10000 }

    output:
    file("${id}.tab")

    // Requires Blast 2.2 separate env
    """
    kleborate -a $assembly --all -o ${id}.tab
    """

  }

}
