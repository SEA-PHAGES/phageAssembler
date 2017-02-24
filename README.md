# phageAssembler
A python script to do some of the repetitive steps of phage genome assembly and quality control.

## What it Does

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

This script uses Newbler (aka GS De Novo Assembler), local blast from NCBI, and AceUtil to perform various steps of the assembly/QC process.  All of these are already installed on the 2017 Science Education Alliance Virtual Machine, maintained and distributed by the University of Pittsburgh and Howard Hughes Medical Institute.  If you are in the SEA-PHAGES program, you have access to download the 2017 (or later) SEA VM [here](http://seaphages.org/software/virtualmachine/).  If you are not a member of the SEA-PHAGES program, but are interested in using this script for academic purposes, [contact us](http://seaphages.org/contact/) to request help.

In addition to what's already included in the SEA VM, you'll need biopython.  To get it, open the seafaculty account on the SEA VM, open a Terminal, and type or paste the following commands.  Enter the seafaculty password if prompted, and select default options if prompted by pressing Enter.

```
sudo apt-get install python-dev
sudo pip install biopython
```

### Installing

Within the 2017 SEA Ubtuntu Virtual Machine, open a terminal.  You can do this by locating the terminal icon from the Launch Bar, or by clicking on the top item in the Launch Bar and searching for "Terminal".

In the terminal type/paste the following command.

```
git clone https://github.com/Danos2000/phageAssembler
```

The code will be copied to a new folder called phageAssembler in your home directory.

Though optional, if you'd like to easily launch the program from any directory in the future, we recommend running the following command in your terminal.

```
sudo cp ~/phageAssembler/phageAssembler.py /usr/local/bin/
```

To check if everything worked, try running the program by typing the following command.

```
phageAssembler.py
```

If it works, you should get a "usage" message that ends in an error for having too few arguments.  All good, you just gave it no input!

If that doesn't work, perhaps the copying didn't work.  Try the following command instead.

```
python2.7 ~/phageAssembler/phageAssembler.py
```

Again, it should result in an error for too few arguments.

### Usage



### Issues

Be warned that this is just a hacked-together script that may or may not work.  Some brief admonitions:

* Do not try to run this script from within a shared folder on your VM.  It will likely fail as it won't have proper permissions to write in the directory.
* Newbler assembly should take only a few minutes, depending on your system, for 50000 reads.  If you use a very large number of reads, you may grind to a halt.

## Authors

* **Dan Russell** - [PhagesDB](http://phagesdb.org/), [seaphages.org](http://seaphages.org/), and the [University of Pittsburgh](http://www.biology.pitt.edu/person/daniel-russell)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks to Becky Garlena for testing and troubleshooting.
* Thanks to Steve Cresawn, Charlie Bowman, and Chris Shaffer for help developing the SEA VM.

