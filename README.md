# MetaCompare

MetaCompare is a computational pipeline for prioritizing resistome risk by estimating the potential for ARGs to be disseminated into human pathogens from a given environmental sample based on metagenomic sequencing data.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine or server.

### System Requirements (*tested on linux Ubuntu 14.04*)

* git installed
* Python3 with pandas package installed
* Blast 2.2.8 or higher version installed

### Installing

**Step 1:** Change the current working directory to the location where you want the cloned `MetaCompare` directory to be made.
**Step 2:** Clone the repository using git command
```
~$ git clone https://github.com/minoh0201/MetaCompare
```

**Step 3:** make directory `BlastDB` and change woring directory to it

```
~$ mkdir BlastDB
~$ cd BlastDB
```

**Step 4:** download the compressed Blast Database file from the web server.

```
~/BlastDB$ wget http://bench.cs.vt.edu/ftp/data/metacomp/BlastDB.tar
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
