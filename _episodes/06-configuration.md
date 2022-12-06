---
title: "Nextflow configuration"
teaching: 30
exercises: 15
questions:
- "What is the difference between the workflow implementation and the workflow configuration?"
- "How do I configure a Nextflow workflow?"
- "How do I assign different resources to different processes?"
- "How do I separate and provide configuration for different computational systems?"
- "How do I change configuration settings from the default settings provided by the workflow?"
objectives:
- "Understand the difference between workflow implementation and configuration."
- "Understand the difference between configuring Nextflow and a Nextflow script."
- "Create a Nextflow configuration file."
- "Understand what a configuration scope is."
- "Be able to assign resources to a process."
- "Be able to refine configuration settings using process selectors."
- "Be able to group configurations into profiles for use with different computer infrastructures."
- "Be able to override existing settings."
- "Be able to inspect configuration settings before running a workflow."
keypoints:
- "Nextflow configuration can be managed using a Nextflow configuration file."
- "Nextflow configuration files are plain text files containing a set of properties."
- "You can define process-specific settings, such as cpus and memory, within the `process` scope."
- "You can assign different resources to different processes using the process selectors `withName` or `withLabel`."
- "You can define a profile for different configurations using the `profiles` scope. These profiles can be selected when launching a pipeline execution by using the `-profile` command-line option"
- "Nextflow configuration settings are evaluated in the order they are read-in."
- "Workflow configuration settings can be inspected using `nextflow config <script> [options]`."
---

## Nextflow configuration

A key Nextflow feature is the ability to decouple the workflow implementation, which describes
the flow of data and operations to perform on that data, from the configuration settings required by the underlying execution platform. This enables the workflow to be portable, allowing it to run on different computational platforms such as an institutional HPC or cloud infrastructure, without needing to modify the workflow implementation.

We have seen earlier that it is possible to provide a `process` with
directives. These directives are process specific configuration settings.
Similarly, we have also provided parameters to our workflow which
are parameter configuration settings. These configuration settings
can be separated from the workflow implementation, into a
configuration file.

## Configuration files

Settings in a configuration file are sets of name-value pairs
(`name = value`). The `name` is a specific property to set,
while the `value` can be anything you can assign to a variable (see
[nextflow scripting](02-nextflow_scripting)), for example, strings,
booleans, or other variables.
It is also possible to access any variable defined in the
host environment such as `$PATH`, `$HOME`, `$PWD`, etc.

~~~
// nextflow.config
my_home_dir = "$HOME"
~~~
{: .language-groovy}

> ## Accessing variables in your configuration file
>
> Generally, variables and functions defined in a
> configuration file are not accessible from the
> workflow script. Only variables defined using the
> `params` scope and the `env` scope (without `env` prefix) can
> be accessed from the workflow script.
>
> ~~~
> workflow {
>     MY_PROCESS( params.input )
> }
> ~~~
> {: .language-groovy}
{: .callout}

Settings are also partitioned into
scopes, which govern the behaviour of different elements of the
workflow. For example, workflow parameters are governed from the
`params` scope, while process directives are governed from the `process` scope. A full list of the available scopes can be found in the
[documentation](https://www.nextflow.io/docs/latest/config.html#config-scopes). It is also possible to define your own scope.

Configuration settings for a workflow are often stored in the file
`nextflow.config` which is in the same directory as the workflow script.
Configuration can be written in either of two ways. The first is using
dot notation, and the second is using brace notation. Both forms
of notation can be used in the same configuration file.

An example of dot notation:
~~~
params.input = ''             // The workflow parameter "input" is assigned an empty string to use as a default value
params.outdir = './results'   // The workflow parameter "outdir" is assigned the value './results' to use by default.
~~~
{: .language-groovy }

An example of brace notation:
~~~
params {
    input  = ''
    outdir = './results'
}
~~~
{: .language-groovy }

Configuration files can also be separated into multiple files and
included into another using the `includeConfig` statement.

~~~
// nextflow.config
params {
    input  = ''
    outdir = './results'
}

includeConfig 'system_resources.config'
~~~
{: .language-groovy}
~~~
// system_resources.config
process {
    cpus = 1    // default cpu usage
    time = '1h' // default time limit
}
~~~
{: .language-groovy}

## How configuration files are combined

Configuration settings can be spread across several files. This also
allows settings to be overridden by other configuration files. The
priority of a setting is determined by the following order,
ranked from highest to lowest.

1. Parameters specified on the command line (`--param_name value`).
1. Parameters provided using the `-params-file` option.
1. Config file specified using the `-c` my_config option.
1. The config file named `nextflow.config` in the current directory.
1. The config file named `nextflow.config` in the workflow project directory (`$projectDir`: the directory where the script to be run is located).
1. The config file `$HOME/.nextflow/config`.
1. Values defined within the workflow script itself (e.g., `main.nf`).

If configuration is provided by more than one of these methods,
configuration is merged giving higher priority to configuration
provided higher in the list.

Existing configuration can be completely ignored by using `-C <custom.config>` to use only configuration provided in the `custom.config` file.

> ## Configuring Nextflow vs Configuring a Nextflow workflow
>
> Parameters starting with a single dash `-` (e.g., `-c my_config.config`) are configuration
> options for `nextflow`, while parameters starting with a double
> dash `--` (e.g., `--outdir`) are workflow parameters defined in the `params` scope.
>
> The majority of Nextflow configuration settings must be provided
> on the command-line, however a handful of settings can also
> be provided within a configuration file, such as
> `workdir = '/path/to/work/dir'` (`-w /path/to/work/dir`) or
> `resume = true` (`-resume`), and do not
> belong to a configuration scope.
{: .callout}

> ## Configuring Nextflow vs Configuring a Nextflow workflow
>
> Parameters starting with a single dash `-` (e.g., `-c my_config.config`) are configuration
> options for `nextflow`, while parameters starting with a double
> dash `--` (e.g., `--outdir`) are workflow parameters defined in the `params` scope.
>
> The majority of Nextflow configuration settings must be provided
> on the command-line, however a handful of settings can also
> be provided within a configuration file, such as
> `workdir = '/path/to/work/dir'` (`-w /path/to/work/dir`) or
> `resume = true` (`-resume`), and do not
> belong to a configuration scope.
{: .callout}

> ## Determine script output
> Determine the outcome of the following script executions.
> Given the script `print_message.nf`:
> ~~~
> nextflow.enable.dsl = 2
>
> params.message = 'hello'
>
> workflow {
>     PRINT_MESSAGE(params.message)
> }
>
> process PRINT_MESSAGE {
>     echo true
>
>     input:
>     val my_message
>
>     script:
>     """
>     echo $my_message
>     """
> }
> ~~~
> {: .language-groovy}
> and configuration (`print_message.config`):
> ~~~
> params.message = 'Are you tired?'
> ~~~
> What is the outcome of the following commands?
> 1. `nextflow run print_message.nf`
> 1. `nextflow run print_message.nf --message '¿Que tal?'`
> 1. `nextflow run print_message.nf -c print_message.config`
> 1. `nextflow run print_message.nf -c pring_message.config --message '¿Que tal?'`
>
> > ## Solution
> >
> > 1. 'hello' - Workflow script uses the value in `print_message.nf`
> > 1. '¿Que tal?' - The command-line parameter overrides the script setting.
> > 1. 'Are you tired?' - The configuration overrides the script setting
> > 1. '¿Que tal?' - The command-line parameter overrides both the script and configuration settings.
> {: .solution}
{: .challenge}


## Configuring process behaviour

Earlier we saw that `process` directives allow the specification of
settings for the task execution such as `cpus`, `memory`, `conda`
and other resources in the pipeline script. This is useful when
prototyping a small workflow script, however this ties the configuration
to the workflow, making it less portable. A good practice is to
separate the process configuration settings into another file.

The `process` configuration scope allows the setting of any process directives in the Nextflow configuration file.

For example:

~~~
// nextflow.config
process {
    cpus = 2
    memory = 8.GB
    time = '1 hour'
    publishDir = [ path: params.outdir, mode: 'copy' ]
}
~~~
{: .language-groovy }

> ## Unit values
>
> Memory and time duration units can be specified either using a string
> based notation in which the digit(s) and the unit can be separated by
> a space character, or
> by using the numeric notation in which the digit(s) and the unit are
> separated by a dot character and not enclosed by quote characters.
>
>  | String syntax   | Numeric syntax | Value                 |
>  |-----------------|----------------|-----------------------|
>  | '10 KB'         | 10.KB          | 10240 bytes           |
>  | '500 MB'        | 500.MB         | 524288000 bytes       |
>  | '1 min'         | 1.min          | 60 seconds            |
>  | '1 hour 25 sec' | -              | 1 hour and 25 seconds |
>
{: .callout}

These settings are applied to all processes in the workflow. A
process selector can be used to apply the configuration to a
specific process or group of processes.

### Process selectors

The resources for a specific process can be defined using `withName:`
followed by the process name ( either the simple name e.g., `'FASTQC'`,
or the fully qualified name e.g., `'NF_RNASEQ:RNA_SEQ:SAMTOOLS_SORT'`),
and the directives within curly braces.
For example, we can specify different `cpus` and `memory` resources
for the processes `INDEX` and `FASTQC` as follows:

~~~
// process_resources.config
process {
    withName: INDEX {
        cpus = 4
        memory = 8.GB
    }
    withName: FASTQC {
        cpus = 2
        memory = 4.GB
    }
}
~~~
{: .language-groovy }

When a workflow has many processes, it is inconvenient to specify
directives for all processes individually, especially if directives
are repeated for groups of processes. A helpful strategy is to annotate
the processes using the `label` directive (processes can have multiple
labels). The `withLabel` selector then allows the configuration of all
processes annotated with a specific label, as shown below:

~~~
// configuration_process_labels.nf
nextflow.enable.dsl=2

process P1 {

    label "big_mem"

    script:
    """
    echo P1: Using $task.cpus cpus and $task.memory memory.
    """
}

process P2 {

    label "big_mem"

    script:
    """
    echo P2: Using $task.cpus cpus and $task.memory memory.
    """
}

workflow {

    P1()
    P2()

}
~~~
{: .language-groovy}

~~~
// configuration_process-labels.config
process {
    withLabel: big_mem {
        cpus = 16
        memory = 64.GB
    }
}
~~~
{: .language-groovy}

Another strategy is to use process selector expressions. Both
`withName:` and `withLabel:` allow the use of regular expressions
to apply the same configuration to all processes matching a pattern.
Regular expressions must be quoted, unlike simple process names
or labels.

- The `|` matches either-or, e.g., `withName: 'INDEX|FASTQC'`
applies the configuration to any process matching the name `INDEX`
or `FASTQC`.
- The `!` inverts a selector, e.g., `withLabel: '!small_mem'` applies
the configuration to any process without the `small_mem` label.
- The `.*` matches any number of characters, e.g.,
`withName: 'NF_RNASEQ:RNA_SEQ:BAM_SORT:.*'` matches all processes
of the workflow `NF_RNASEQ:RNA_SEQ:BAM_SORT`.

A regular expression cheat-sheet can be found
[here](https://www.jrebel.com/system/files/regular-expressions-cheat-sheet.pdf) if you would like to write more expressive expressions.

#### Dynamic expressions

A common scenario is that configuration settings may depend on the
data being processed. Such settings can be dynamically expressed
using a closure. For example, we can specify the `memory` required
as a multiple of the number of `cpus`. Similarly, we can publish
results to a subfolder based on the sample name.

~~~
process FASTQC {

    input:
    tuple val(sample), path(reads)

    script:
    """
    fastqc -t $task.cpus $reads
    """
}
~~~
{: .language-groovy }

~~~
// nextflow.config
process {
    withName: FASTQC {
        cpus = 2
        memory = { 2.GB * task.cpus }
        publishDir = { "fastqc/$sample" }
    }
}
~~~
{: .language-groovy }

## Configuring execution platforms

Nextflow supports a wide range of execution platforms, from
running locally, to running on HPC clusters or cloud infrastructures.
See https://www.nextflow.io/docs/latest/executor.html for the
full list of supported executors.

![nf-executors](https://seqera.io/training/img/nf-executors.png)

The default executor configuration is defined within the `executor`
scope (https://www.nextflow.io/docs/latest/config.html#scope-executor).
For example, in the config below we specify the executor as
Sun Grid Engine, `sge` and the number of tasks the executor will
handle in a parallel manner (`queueSize`) to 10.

~~~
// nextflow.config
executor {
    name = 'sge'
    queueSize = 10
}
~~~
{: .language-groovy }

The `process.executor` directive allows you to override
the executor to be used by a specific process. This can be
useful, for example, when there are short running tasks
that can be run locally, and are unsuitable for submission
to HPC executors (check for guidelines on best practice use
of your execution system). Other process directives such as
`process.clusterOptions`, `process.queue`, and `process.machineType`
can be also be used to further configure processes depending
on the executor used.  

~~~
//nextflow.config
executor {
    name = 'sge'
    queueSize = 10
}
process {
    withLabel: 'short' {
        executor = 'local'
    }
}
~~~
{: .language-groovy }

## Configuring software requirements

An important feature of Nextflow is the ability to manage
software using different technologies. It supports the Conda package
management system, and container engines such as Docker, Singularity,
Podman, Charliecloud, and Shifter. These technologies
allow one to package tools and their dependencies into a software environment
such that the tools will always work as long as the environment can be loaded.
This facilitates portable and reproducible workflows.
Software environment specification is managed from the `process` scope,
allowing the use of process selectors to manage which processes
load which software environment. Each technology also has its own
scope to provide further technology specific configuration settings.

### Software configuration using Docker

Docker is a container technology. Container images are
lightweight, standalone, executable package of software
that includes everything needed to run an application:
code, runtime, system tools, system libraries and settings.
Containerized software is intended to run the same regardless
of the underlying infrastructure, unlike other package management
technologies which are operating system dependant (See
the [published article on Nextflow](https://doi.org/10.1038/nbt.3820)).
For each container image used, Nextflow uses Docker to spawn
an independent and isolated container instance for each process task.

To use Docker, we must provide a container image path using the
`process.container` directive, and also enable docker in the docker
scope, `docker.enabled = true`. A container image path takes the form
`(protocol://)registry/repository/image:version--build`.
By default, Docker containers
run software using a privileged user. This can cause issues,
and so it is also a good idea to supply your user and group
via the `docker.runOptions`.

~~~
process.container = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
~~~
{: .language-groovy }

### Software configuration using Singularity

Singularity is another container technology, commonly used on
HPC clusters. It is different to Docker in several ways. The
primary differences are that processes are run as the user,
and certain directories are automatically "mounted" (made available)
in the container instance. Singularity also supports building
Singularity images from Docker images, allowing Docker image paths
to be used as values for `process.container`.

Singularity is enabled in a similar manner to Docker.
A container image path must be provided using `process.container` and
singularity enabled using `singularity.enabled = true`.

~~~
process.container = 'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0'
singularity.enabled = true
~~~
{: .language-groovy }

> ## Container protocols
>
> The following protocols are supported:
> - `docker://``: download the container image from the Docker Hub and convert it to the Singularity format (default).
> - `library://``: download the container image from the Singularity Library service.
> - `shub://``: download the container image from the Singularity Hub.
> - `docker-daemon://`: pull the container image from a local Docker installation and convert it to a Singularity image file.
> - `https://`: download the singularity image from the given URL.
> - `file://`: use a singularity image on local computer storage.
{: .callout}

## Configuration profiles

One of the most powerful features of Nextflow configuration is to
predefine multiple configurations or `profiles` for different
execution platforms. This allows a group of predefined settings to
be called with a short invocation, `-profile <profile name>`.

Configuration profiles are defined in the `profiles` scope,
which group the attributes that belong to the same profile
using a common prefix.

~~~
//Example1: configuration_profiles.config
profiles {

    standard {
        params.genome = '/local/path/ref.fasta'
        process.executor = 'local'
    }

    cluster {
        params.genome = '/data/stared/ref.fasta'
        process.executor = 'sge'
        process.queue = 'long'
        process.memory = '10GB'
        process.conda = '/some/path/env.yml'
    }

    cloud {
        params.genome = '/data/stared/ref.fasta'
        process.executor = 'awsbatch'
        process.container = 'cbcrg/imagex'
        docker.enabled = true
    }

}
~~~
{: .language-groovy }

<b>Take note</b>: The parameters defined in the nextflow.config file overrides the default parameters inside your nextflow script! However, the command-line parameters even override the nextflow.config file!
By default, we are running local executor but you can configure your processes to use other executors (slurm, sge, awsbatch, google-lifesciences).
Nextflow requires a shared storage to exchange data between tasks.

~~~
// This example defines execution profiles for different environments
// Example2: configuration_profiles.config
profiles {
  slurm {
    process.container = 'nextflow/rnaseq-nf:latest'  // you can use a local container, '/path/to/image.sif'
    singularity.enabled = true
    process.executor = 'slurm'
    process.queue = 'compute'
    process.time = '8h'
    clusterOptions = '-q dev'
    module = ['singularity']
    process.memory = '8 GB'
    process.cpus = 8
  }
  standard {
   singularity.enabled = true
   process.container = 'nextflow/rnaseq-nf:latest' // you can use a local container, '/path/to/image.sif'
   process.executor = 'local'
 }
  awsbatch {
    process.executor = 'awsbatch'
    process.queue = 'nextflow-ci'
    process.container = 'nextflow/rnaseq-nf:latest'
    workDir = 's3://nextflow-ci/work/'
    aws.region = 'eu-west-1'
    aws.batch.cliPath = '/home/ec2-user/miniconda/aws'
    params.reads = 's3://rnaseq-nf/data/ggal/lung_{1,2}.fq'
    params.transcriptome = 's3://rnaseq-nf/data/ggal/transcript.fa'
  }
  's3-data'{
    process.container = 'nextflow/rnaseq-nf:latest'
    params.reads = 's3://rnaseq-nf/data/ggal/lung_{1,2}.fq'
    params.transcriptome = 's3://rnaseq-nf/data/ggal/transcript.fa'    
  }
  gls {
    process.executor = 'google-lifesciences'
    process.container = 'nextflow/rnaseq-nf:latest'
    workDir = 'gs://nfdemo/scratch' // <- replace with your own bucket!
    google.region = 'europe-west2'
    params.transcriptome = 'gs://rnaseq-nf/data/ggal/transcript.fa'
    params.reads = 'gs://rnaseq-nf/data/ggal/gut_{1,2}.fq'
    params.transcriptome = 'gs://rnaseq-nf/data/multiqc'
  }
  conda {
    process.conda = "$baseDir/conda.yml"
  }
}

aws {
    accessKey = '<YOUR S3 ACCESS KEY>'
    secretKey = '<YOUR S3 SECRET KEY>'
    region = '<AWS REGION IDENTIFIER>'
}
google {
    project = 'YOUR PROJECT ID HERE'
    location = 'us-central1'
    accessKey = 'YOUR API KEY HERE'
}

~~~
{: .language-groovy }

To run your Nexflow script with awsbatch, you can use:

~~~
$ nextflow run main.nf -profile awsbatch
~~~
{: .language-bash}

This configuration defines different profiles: `standard`,
`slurm`, `awsbatch`, `s3-data`, `gls` and `conda` that set different process configuration
strategies depending on the target execution platform. By
convention the standard profile is implicitly used when no
other profile is specified by the user. To enable a specific
profile use `-profile` option followed by the profile name:

~~~
nextflow run <your script> -profile standard
~~~
{: .language-bash}

~~~
nextflow run wc.nf -profile standard
~~~
{: .language-bash}

~~~
N E X T F L O W  ~  version 22.10.1
Launching `wc.nf` [nasty_hilbert] DSL2 - revision: ed06b3439b
executor >  local (1)
[f3/560c2f] process > NUM_LINES (1) [100%] 1 of 1 ✔
ref1.fa
Number of lines: 2852 and Number of cpus: 2
my_script -m 2 GB -n 2 -t 1h
~~~
{: .output}

> ## Configuration order
>
> Settings from profiles will override general settings
> in the configuration file. However, it is also important
> to remember that configuration is evaluated in the order it
> is read in. For example, in the following example, the `publishDir`
> directive will always take the value 'results' even when the
> profile `hpc` is used. This is because the setting is evaluated
> before Nextflow knows about the `hpc` profile. If the `publishDir`
> directive is moved to after the `profiles` scope, then `publishDir`
> will use the correct value of `params.results`.
>
> ~~~
> params.results = 'results'
> process.publishDir = params.results
> profiles {
>     hpc {
>         params.results = '/long/term/storage/results'
>     }
> }
> ~~~
> {: .language-groovy}
{: .callout}

## Inspecting the Nextflow configuration

You can use the command `nextflow config` to print the resolved
configuration of a workflow. This allows you to see what settings
Nextflow will use to run a workflow.

~~~
$ nextflow config wc.nf -profile slurm
~~~
{: .language-bash}

~~~
singularity {
   enabled = true
}

executor = 'slurm'
queue = 'compute'
time = '8h'
clusterOptions = '-q dev'
module = ['singularity']
memory = '8 GB'
cpus = 8

~~~
{: .output}



## Additional Nextflow things

<b>Nextflow Logging:</b>
~~~
$ nextflow log
~~~
{: .language-bash}

~~~
TIMESTAMP          	DURATION	RUN NAME            	STATUS	REVISION ID	SESSION ID                          	COMMAND             
2022-11-26 03:39:58	13s     	extravagant_gates   	OK    	75c1898388 	36b63809-36cc-4b1c-acd8-5ce5becf13fb	nextflow run main.nf
2022-12-05 16:25:29	3s      	pedantic_wescoff    	ERR   	e99fe11539 	a0876c71-eec9-43ee-ac57-f924404a19e0	nextflow run wc.nf  
2022-12-05 16:27:49	2.6s    	dreamy_roentgen     	OK    	e99fe11539 	0bef4823-caeb-4642-a012-9563cacd938d	nextflow run wc.nf  
2022-12-05 16:51:27	2.8s    	hungry_kowalevski   	ERR   	7899ed9c6a 	631044c7-1ac6-41dc-a690-635d46ee27da	nextflow run wc.nf
~~~
{: .output}

~~~
$ nextflow log extravagant_gates -F 'process == /multiqc/'
~~~
{: .language-bash}

~~~
/home/hpc_user/Desktop/nfdemo/work/83/9ec8858202438d3b7097d8c26912c4
~~~
{: .output}

<b>Running Nextflow Pipelines (GitHub):</b>
~~~
$ nextflow run https://github.com/nextflow-io/hello
~~~
{: .language-bash}

~~~
N E X T F L O W  ~  version 22.10.1
Launching `https://github.com/nextflow-io/hello` [clever_euclid] DSL2 - revision: 4eab81bd42 [master]
executor >  local (4)
[99/92026c] process > sayHello (2) [100%] 4 of 4 ✔
Hola world!

Bonjour world!

Hello world!

Ciao world!
~~~
{: .output}

~~~
$ nextflow run nextflow-io/hello -profile standard
~~~
{: .language-bash}

~~~
N E X T F L O W  ~  version 22.10.1
Launching `https://github.com/nextflow-io/hello` [tender_khorana] DSL2 - revision: 4eab81bd42 [master]
executor >  local (4)
[21/257184] process > sayHello (3) [100%] 4 of 4 ✔
Bonjour world!

Ciao world!

Hola world!

Hello world!
~~~
{: .output}

<b>Task Execution Caching:</b>
~~~
$ nextflow run main.nf
~~~
{: .language-bash}

~~~
N E X T F L O W  ~  version 22.10.1
Launching `main.nf` [cheeky_raman] DSL2 - revision: 75c1898388
RNASEQ NEXTFLOW PIPELINE uses the ffg bioinformatics tools: Salmon, FastQC, MultiQC
transcriptome : /home/hpc_user/Desktop/hpc_user/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa
reads         : /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_{1,2}.fq
outdir        : results

executor >  local (4)
[57/3676c0] process > rnaseq_sub:index (ggal_1_48850000_49020000) [100%] 1 of 1 ✔
[ce/017ea6] process > rnaseq_sub:fastqc (FASTQC on ggal_gut)      [100%] 1 of 1 ✔
[87/78e582] process > rnaseq_sub:quant (1)                        [100%] 1 of 1 ✔
[d0/a66866] process > multiqc                                     [100%] 1 of 1 ✔
[ggal_gut, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/hpc_user/Desktop/hpc_user/data/ggal/ggal_gut_2.fq]
[ggal_gut, [/home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/hpc_user/Desktop/hpc_user/data/ggal/ggal_gut_2.fq]]
~~~
{: .output}

~~~
$ nextflow run main.nf --reads 'data/ggal/ggal_*_{1,2}.fq' -resume
~~~
{: .language-bash}

~~~
N E X T F L O W  ~  version 22.10.1
Launching `main.nf` [disturbed_panini] DSL2 - revision: 75c1898388
RNASEQ NEXTFLOW PIPELINE uses the ffg bioinformatics tools: Salmon, FastQC, MultiQC
transcriptome : /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa
reads         : data/ggal/ggal_*_{1,2}.fq
outdir        : results

executor >  local (3)
[87/4fba91] process > rnaseq_sub:index (ggal_1_48850000_49020000) [100%] 1 of 1, cached: 1 ✔
[9a/92a372] process > rnaseq_sub:fastqc (FASTQC on ggal_liver)    [100%] 2 of 2, cached: 1 ✔
[32/c6022a] process > rnaseq_sub:quant (2)                        [100%] 2 of 2, cached: 1 ✔
[aa/57ea75] process > multiqc                                     [100%] 1 of 1 ✔
[ggal_gut, [/home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]]
[ggal_gut, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_1.fq, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_gut_2.fq]
[ggal_liver, [/home/hpc_user/Desktop/nfdemo/data/ggal/ggal_liver_1.fq, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_liver_2.fq]]
[ggal_liver, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_liver_1.fq, /home/hpc_user/Desktop/nfdemo/data/ggal/ggal_liver_2.fq]
~~~
{: .output}

{% include links.md %}
