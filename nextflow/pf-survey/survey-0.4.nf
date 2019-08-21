log.info """
=======================================================
                     PathFinder
                Survey Pipeline v0.4
               (QC, Assembly, Typing)
=======================================================

Directory configuration:

Output directory:               $params.outdir
Fastq input:                    $params.fastq
Resource directory:             $params.resources

Process configuration:

leading         =    $params.leading
trailing        =    $params.trailing
min_len         =    $params.min_len
sliding_window  =    $params.sliding_window
adapters        =    $params.adapters
kraken_db       =    $params.taxdb
assembler       =    $params.assembler
depth           =    $params.depth
reference       =    $params.reference

mash            =    $params.mash
snippy          =    $params.snippy
mykrobe         =    $params.mykrobe
kleborate       =    $params.kleborate

=======================================================
=======================================================
"""


fastq = Channel
		.fromFilePairs("${params.fastq}", flat: true)
		.ifEmpty { exit 1, "Input read files could not be found." }

if (params.genus) {
    use_genus = "--usegenus"
} else {
    use_genus = ""
}

process Trimmomatic {

  label "trimmomatic"

	tag { id }

	input:
	set id, file(forward), file(reverse) from fastq

	output:
	set id, file("${id}_1.fq.gz"), file("${id}_2.fq.gz") into kraken, snippy, shovill, mykrobe

	"""
	trimmomatic PE ${forward} ${reverse} -threads $task.cpus $params.phred -baseout ${id}.fq.gz \
	ILLUMINACLIP:${params.resources}/${params.adapters}:2:30:10:3:TRUE LEADING:${params.leading} TRAILING:${params.trailing} \
	SLIDINGWINDOW:${params.sliding_window} MINLEN:${params.min_len}

	mv ${id}_1P.fq.gz ${id}_1.fq.gz
	mv ${id}_2P.fq.gz ${id}_2.fq.gz
	"""

}

process Shovill {

    label "shovill"

    tag { id }
    publishDir "${params.outdir}/skesa", mode: "copy"

    input:
    set id, file(forward), file(reverse) from shovill.filter{ it[1].size()>10000 & it[2].size()>10000 } // Failed / empty file can be up to 10kb

    output:
    set id, file("${id}.fasta") into mash, abricate, mlst, prokka, kleborate, sccion

    """
    shovill --R1 ${forward} --R2 ${reverse} --cpus $task.cpus --ram $task.memory \
    --depth $params.depth --assembler $params.assembler --outdir $id --force
    mv ${id}/contigs.fa ${id}.fasta
    """

}


process Kraken2 {

    label "kraken2"

    tag { id }
    publishDir "${params.outdir}/kraken", mode: "copy"

    input:
    set id, file(forward), file(reverse) from kraken.filter{ it[1].size()>10000 & it[2].size()>10000 }

    output:
    set id, file("${id}.report")

    """

    kraken2 --db ${params.resources}/${params.minikraken} --threads $task.cpus --output "-" \
    --fastq-input --gzip-compressed --paired --report ${id}.report --use-names ${forward} ${reverse}
    """

}

if (params.snippy){

    process Snippy {

      label "snippy"
      tag { id }
      publishDir "${params.outdir}/snippy", mode: "copy"

      input:
      set id, file(forward), file(reverse) from snippy

      output:
      file("$id/${id}.tab") into core_alignment

      """
      snippy --cpus $task.cpus --ram $task.memory --outdir $id --prefix $id \
      --reference $params.resources/$params.reference --R1 $forward --R2 $reverse
      """

    }

}

if (params.mykrobe) {
    process MykrobePredictor {

        label "mykrobe"

        tag { id }
        publishDir "${params.outdir}/mykrobe", mode: "copy"

        input:
        set id, file(forward), file(reverse) from mykrobe.filter{ it[1].size()>10000 & it[2].size()>10000 }

        output:
        set id, file("${id}.json")

        """
        mykrobe predict --format json --threads $task.cpus ${id} $params.mykrobe_other --seq $forward $reverse > ${id}.json
        """

    }
}

process Abricate {

  label "abricate"

	// Assembly-based typing with T. Seemann's Abricate

	tag { id }

	publishDir "${params.outdir}/abricate", mode: "copy"

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

  label "mlst"

	// Assembly-based MLST

	tag { id }
	publishDir "${params.outdir}/mlst", mode: "copy"

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

      label "mash"

      // MinHash distance against global reference:

      tag { id }
      publishDir "${params.outdir}/mash", mode: "copy"

      input:
      set id, file(assembly) from mash.filter{ it[1].size()>10000 }

      output:
      file("${id}.mash.tab")

      """
      mash dist $assembly $params.resources/$params.reference > ${id}.mash.tab
      """
    }
}

process Prokka {

  label "prokka"

  // MinHash distance against global reference:

  tag { id }
  publishDir "${params.outdir}/prokka", mode: "copy"

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

    label "kleborate"

    // MinHash distance against global reference:

    tag { id }
    publishDir "${params.outdir}/kleborate", mode: "copy"

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

if (params.sccion) {

  process SCCion {

    label "sccion"

    // SCCion typing de novo assembly:

    tag { id }
    publishDir "${params.outdir}/sccion", mode: "copy"

    input:
    set id, file(assembly) from sccion.filter{ it[1].size()>10000 }

    output:
    file("${id}.tab")

    """
    sccion type -a $assembly --resistance --mge > ${id}.tab
    """

  }

}
