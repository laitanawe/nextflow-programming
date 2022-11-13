---
title: "Nextflow DSL Workflow Definition"
teaching: 50
exercises: 10
questions:
- "What are the main differences between a Nextflow DSL2 and DSL1 process?"
- "What are the main components of a Nextflow DSL2 pipeline?"
- "How do I run a Nextflow pipeline?"
objectives:
- "Explain the benefits of using Nextflow DSL2 to develop pipelines."
- "Explain the components of a Nextflow DSL2 process, pipeline."
- "Run a Nextflow DSL2 pipeline."
keypoints:
- "A DSL2 workflow is a sequence of tasks that process a set of data."
- "In Nextflow DSL2 processes, the process definition is no longer tied to specific channels and so, it can be used independently. This enables reuse of DSL2 components. From some module, we can include DSL2 processes in a workflow."
- "You run Nextflow DSL2 modules and pipelines using the `nextflow run` command."
---

## Our example script

We are now going to look at a sample Nextflow DSL2 workflow that performs some RNA-seq tasks.

For the Nextflow course demo, you need to download some files to follow the lesson.
1. Download <a href="https://laitanawe.github.io/nextflow-novice/data/nextflow-nov-lesson.tar.gz">nextflow-nov-lesson.tar.gz</a>
2. Make a directory called "nfdemo" on the Desktop
3. Move the nextflow-nov-lesson.tar.gz inside the nfdemo directory and cd into the nfdemo directory.
4. Unzip/extract `nextflow-nov-lesson.tar.gz` by typing `tar -xzvf nextflow-nov-lesson.tar.gz`
You should end up certain files within the folder **`nfdemo/data/ggal`** on your Desktop.

Open the file `main.nf` in the script directory with your favourite text editor.

This is a Nextflow DSL2 pipeline performing the following steps:
1. Indexes a transcriptome file.
1. Performs quality control.
1. Performs quantification.
1. Creates a MultiQC report.

~~~
#!/usr/bin/env nextflow

/* Contributors:
 * - Awe, O.I
 */

nextflow.enable.dsl = 2

include { index } from './modules/index' // './modules/some_module' //Cut and Paste index block into ./modules/index.nf
include { quant } from './modules/quant' // './modules/other_module' //Cut and Paste index block into ./modules/quant.nf
include { fastqc } from './modules/fastqc' // './modules/another_more_module' //Cut and Paste index block into ./modules/fastqc.nf
include { multiqc } from './modules/multiqc' // './modules/one_more_module' //Cut and Paste index block into ./modules/multiqc.nf

//  The default workflow
workflow {

transcriptome = params.transcriptome
read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists: true )

    index( transcriptome )
    fastqc( read_pairs_ch )
    quant( index.out, read_pairs_ch )
    multiqc( fastqc.out.mix( quant.out ).collect() )

}
~~~~
{: .language-groovy}

To run a Nextflow script use the command `nextflow run <script_name>`.

> ## Run a Nextflow  script
> Run the script by entering the following command in your terminal:
>
> ~~~
> $ nextflow run main.nf
> ~~~
> {: .language-bash}
> > ## Solution
> > You should see output similar to the text shown below:
> >
> > ~~~
> > N E X T F L O W  ~  version 20.10.0
> > Launching `main.nf` [sleepy_avogadro] DSL2 - revision: 83e5d597be
> > RNASEQ NEXTFLOW PIPELINE uses the ffg bioinformatics tools: Salmon, FastQC, MultiQC
> > transcriptome : /home/aweo/Desktop/nfdemo/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa
> > reads         : /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_{1,2}.fq
> > outdir        : results
> > executor >  local (4)
> > [3b/a13c45] process > index (ggal_1_48850000_49020000) [100%] 1 of 1 ✔
> > [ff/f60afe] process > fastqc (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
> > [34/84f09d] process > quant (1)                        [100%] 1 of 1 ✔
> > [aa/25f8e4] process > multiqc                          [100%] 1 of 1 ✔
> > [ggal_gut, [/home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]]
> > [ggal_gut, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/aweo/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]
> > ~~~
> > {: .output}
> >
> > 1. The first line shows the Nextflow version number.
> > 1. The second line shows the run name `fervent_babbage` (adjective and scientist name) and revision id `c54a707593`.
> > 1. The third line tells you the process has been executed locally (`executor >  local`).
> > 1. The next line shows the log info followed by the process ids e.g. `3b/a13c45`, process name, number of cpus, percentage task completion, and how many instances of the process have been run.
> > 1. The final line is the output of the `view` operator.
> {: .solution}
{: .challenge}


> ## Process identification
> The hexadecimal numbers, like 3b/a13c45, identify the unique process execution.
> These numbers are also the prefix of the directories where each process is executed.
> You can inspect the files produced by changing to the directory `$PWD/work` and
> using these numbers to find the process-specific execution path. We will learn how exactly
> nextflow using *work* directory to execute processes in the following sections.
{: .callout}



{% include links.md %}
