In general to submit a job

~~~shell
qsub -cwd -q std.q -N job_name job_script [ argument ... ]
~~~

but it is recomended to use 

~~~shell
qsub job_script [ argument ... ]
~~~

and add the qsub options directly in the shell script like this:

~~~
#!/bin/bash

#$ -q std.q
#$ -cwd
#$ -N job_name

./your_commands
~~~

The lines that begin with #$ are parsed by the qsub command but are treated as comments by the shell

The useful options are:

- **-q :** used to select the queue, std.q is the standard queue
- **-cwd : ** is one of the most important, this commad makes so the script is executed in the workong directory instead of $HOME
- **-N :**  Name of the job, if not mentioned its the file name of the script
- **-o and -e :** standard output/erro will be appended to the mentioned file 
- **-M :** specifies the email where to send notifications
- **-t m-n:** Trivial parallelization using job arrays

### Examples

~~~
qsub -q high -cwd simple.sh # lounch a job
qdel -u <user_name> # delete all jobs from one user
~~~



